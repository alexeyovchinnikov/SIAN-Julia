struct ODE{P}
    poly_ring::MPolyRing
    x_vars::Array{P, 1}
    y_vars::Array{P, 1}
    u_vars::Array{P, 1}
    parameters::Array{P, 1}
    x_equations::OrderedDict{P, <: Union{P, Generic.Frac{P}}}
    y_equations::OrderedDict{P, <: Union{P, Generic.Frac{P}}}
    
    function ODE{P}(
            x_eqs::OrderedDict{P, <: Union{P, Generic.Frac{P}}}, 
            y_eqs::OrderedDict{P, <: Union{P, Generic.Frac{P}}},    
            inputs::Array{P, 1}
        ) where {P <: MPolyElem{<: FieldElem}}
        # Initialize ODE
        # x_eqs is a dictionary x_i => f_i(x, u, params)
        # y_eqs is a dictionary y_i => g_i(x, u, params)

        num, den = unpack_fraction(collect(values(x_eqs))[1])
        poly_ring = parent(num)
        x_vars = collect(keys(x_eqs))
        y_vars = collect(keys(y_eqs))
        u_vars = inputs
        parameters = filter(v -> (!(v in x_vars) && !(v in u_vars) && !(v in y_vars)), gens(poly_ring))
        new(poly_ring, x_vars, y_vars, u_vars, parameters, x_eqs, y_eqs)
    end
end

#------------------------------------------------------------------------------

function set_parameter_values(ode::ODE{P}, param_values::OrderedDict{P, T}) where {T <: FieldElem, P <: MPolyElem{T}}
    """
    Input:
        - ode, an ODE as above
        - param_values, values for (some of) the parameters as dictionary parameter => value
    Output: new ode with the parameters in param_values plugged with the given numbers
    """
    new_vars = map(var_to_str, [v for v in gens(ode.poly_ring) if !(v in keys(param_values))])
    small_ring, small_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_vars)
    eval_dict = OrderedDict(str_to_var(v, ode.poly_ring) => str_to_var(v, small_ring) for v in new_vars)
    merge!(eval_dict, OrderedDict(p => small_ring(val) for (p, val) in param_values))

    return ODE{P}(
        OrderedDict{P, Union{P, Generic.Frac{P}}}(eval_at_dict(v, eval_dict) => eval_at_dict(f, eval_dict) for (v, f) in ode.x_equations),
         OrderedDict{P, Union{P, Generic.Frac{P}}}(eval_at_dict(v, eval_dict) => eval_at_dict(f, eval_dict) for (v, f) in ode.y_equations),
        [eval_at_dict(u, eval_dict) for u in ode.u_vars]
    )
end

#------------------------------------------------------------------------------

function print_for_SIAN(ode::ODE{P}, outputs::Array{P, 1}) where P <: MPolyElem{<: FieldElem}
    """
    Prints the ODE in the format accepted by SIAN (https://github.com/pogudingleb/SIAN)
    """
    vars_str = OrderedDict(x => var_to_str(x) * "(t)" for x in vcat(ode.x_vars, ode.u_vars))
    merge!(vars_str, OrderedDict(p => var_to_str(p) for p in ode.parameters))
    R_print, vars_print = Nemo.PolynomialRing(base_ring(ode.poly_ring), [vars_str[v] for v in gens(ode.poly_ring)])
    result = ""

    function _lhs_to_str(lhs)
        num, den = unpack_fraction(lhs)
        result = string(evaluate(num, vars_print))
        if den != 1
             result = "($result) / ($(evaluate(den, vars_print)))"
        end
        return result
    end

    for (x, f) in ode.equations
        result = result * "diff(" * var_to_str(x) * "(t), t) = $(_lhs_to_str(f)), \n"
    end
    for (y_ind, g) in enumerate(outputs)
        result = result * "y_var_$y_ind(t) = $(_lhs_to_str(g)), \n"
    end
    return result
end

#------------------------------------------------------------------------------

function macrohelper_extract_vars(equations::Array{Expr, 1})
    funcs, x_vars, all_symb = Array{Any}(undef, 0), Array{Any}(undef, 0), Array{Any}(undef, 0)
    aux_symb = Set([:(+), :(-), :(=), :(*), :(^), :t, :(/), :(//)])
    for eq in equations
        MacroTools.postwalk(
            x -> begin 
                if @capture(x, f_'(t)) 
                    push!(x_vars, f)
                    push!(all_symb, f)
                elseif @capture(x, f_(t))
                    push!(funcs, f)
                elseif (x isa Symbol) && !(x in aux_symb)
                    push!(all_symb, x)
                end
                return x
            end, 
            eq
        )
    end
    io_vars = setdiff(funcs, x_vars)
    all_symb = vcat(x_vars, io_vars, setdiff(all_symb, funcs))
    return x_vars, io_vars, all_symb
end

#------------------------------------------------------------------------------

function macrohelper_clean(ex::Expr)
    ex = MacroTools.postwalk(x -> @capture(x, f_'(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> @capture(x, f_(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> x == :(/) ? :(//) : x, ex)
    return ex
end

#------------------------------------------------------------------------------

macro ODEmodel(ex::Expr...)
    """
    Macros for creating an ODE from a list of equations
    Also injects all variables into the global scope
    """
    equations = [ex...]
    x_vars, io_vars, all_symb = macrohelper_extract_vars(equations)
    
    # creating the polynomial ring
    vars_list = :([$(all_symb...)])
    R = gensym()
    vars_aux = gensym()
    exp_ring = :(($R, $vars_aux) = Nemo.PolynomialRing(Nemo.QQ, map(string, $all_symb)))
    assignments = [:($(all_symb[i]) = $vars_aux[$i]) for i in 1:length(all_symb)]
    
    # preparing equations
    equations = map(macrohelper_clean, equations)
    x_dict = gensym()
    y_dict = gensym()
    y_vars = Array{Any}(undef, 0)
    x_dict_create_expr = :($x_dict = OrderedDict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}())
    y_dict_create_expr = :($y_dict = OrderedDict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}())
    eqs_expr = []
    for eq in equations
        if eq.head != :(=)
            throw("Problem with parsing at $eq") 
        end
        lhs, rhs = eq.args[1:2]
        loc_all_symb = macrohelper_extract_vars([rhs])[3]
        to_insert = undef
        if lhs in x_vars
            to_insert = x_dict
        elseif lhs in io_vars
            to_insert = y_dict
            push!(y_vars, lhs)
        else
            throw("Unknown left-hand side $lhs")
        end
        if isempty(loc_all_symb)
            push!(eqs_expr, :($to_insert[$lhs] = $R($rhs)))
        else
            push!(eqs_expr, :($to_insert[$lhs] = ($rhs)))
        end
    end

    u_vars = setdiff(io_vars, y_vars)
    params = setdiff(all_symb, union(x_vars, y_vars, u_vars))
    print("Summary of the model:\n")
    print("State variables: [", join(map(string, x_vars), ", "), "]\n")
    print("Parameter: [", join(map(string, params), ", "), "]\n")
    print("Inputs: [", join(map(string, u_vars), ", "), "]\n")
    print("Outputs: [", join(map(string, y_vars), ", "), "]\n")
   
    # creating the ode object
    ode_expr = :(ODE{fmpq_mpoly}($x_dict, $y_dict, Array{fmpq_mpoly}([$(u_vars...)])))
    
    result = Expr(
        :block, 
        exp_ring, assignments..., 
        x_dict_create_expr, y_dict_create_expr, eqs_expr..., 
        ode_expr
    )
    return esc(result)
end

#------------------------------------------------------------------------------

function generate_replica(ode::ODE{P}, r::Int) where P <: MPolyElem
    """
    Returns ode_r, and r-fold replica of the original ode.
    States, outputs, and inputs are replicated, parameters are not
    """
    new_varnames = Array{String}(undef, 0)
    for v in vcat(ode.x_vars, ode.y_vars, ode.u_vars)
        append!(new_varnames, [var_to_str(v) * "_r$i" for i in 1:r])
    end
    append!(new_varnames, map(string, ode.parameters))
    new_ring, new_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_varnames)
    new_x_eqs = OrderedDict{P, Union{P, Generic.Frac{P}}}()
    new_y_eqs = OrderedDict{P, Union{P, Generic.Frac{P}}}()
    new_us = Array{P, 1}()
    for i in 1:r
        eval = merge(
            OrderedDict(v => str_to_var(var_to_str(v) * "_r$i", new_ring) for v in vcat(ode.x_vars, ode.y_vars, ode.u_vars)),
            OrderedDict(p => switch_ring(p, new_ring) for p in ode.parameters)
        )
        eval_vec = [eval[v] for v in gens(ode.poly_ring)]
        new_x_eqs = merge(
            new_x_eqs, 
            OrderedDict{P, Union{P, Generic.Frac{P}}}(evaluate(x, eval_vec) => evaluate(f, eval_vec) for (x, f) in ode.x_equations)
        )
        new_y_eqs = merge(
            new_y_eqs,
            OrderedDict{P, Union{P, Generic.Frac{P}}}(evaluate(x, eval_vec) => evaluate(f, eval_vec) for (x, f) in ode.y_equations)
        )
        append!(new_us, [str_to_var(var_to_str(u) * "_r$i", new_ring) for u in ode.u_vars])
    end
    return ODE{P}(new_x_eqs, new_y_eqs, new_us)
end

#------------------------------------------------------------------------------

function _reduce_poly_mod_p(poly::MPolyElem{Nemo.fmpq}, p::Int)
    """
    Reduces a polynomial over Q modulo p
    """
    den = denominator(poly)
    num = change_base_ring(Nemo.ZZ, den * poly)
    if Nemo.GF(p)(den) == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
    end
    return change_base_ring(Nemo.GF(p), num) * (1 // Nemo.GF(p)(den))
end

#------------------------------------------------------------------------------

