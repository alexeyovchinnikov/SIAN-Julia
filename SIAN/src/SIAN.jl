module SIAN
ENV["NEMO_PRINT_BANNER"] = "false"
using Oscar
using LinearAlgebra
using Singular
using GroebnerBasis
using MacroTools
using OrderedCollections

export ODE, @ODEmodel, macroexpand, macrohelper_extract_vars, macrohelper_clean, Nemo, OrderedDict, fmpq_mpoly, Generic, get_parameters, identifiability_ode

struct ODE{P}
    poly_ring::MPolyRing
    x_vars::Array{P,1}
    y_vars::Array{P,1}
    u_vars::Array{P,1}
    parameters::Array{P,1}
    x_equations::OrderedDict{P,<: Union{P,Generic.Frac{P}}}
    y_equations::OrderedDict{P,<: Union{P,Generic.Frac{P}}}

    function ODE{P}(
            x_eqs::OrderedDict{P,<: Union{P,Generic.Frac{P}}},
            y_eqs::OrderedDict{P,<: Union{P,Generic.Frac{P}}},
            inputs::Array{P,1}
        )  where {P <: MPolyElem{<: FieldElem}}
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

    """
        func eval_at_dict(poly::P, d::OrderedDict{P,<: RingElem}) where P <: MPolyElem
        
    Evaluates a polynomial on a dict `var => val` missing values are replaced with zeroes
    """
function eval_at_dict(poly::P, d::OrderedDict{P,<: RingElem}) where P <: MPolyElem
    
    point = [get(d, v, base_ring(parent(poly))(0)) for v in gens(parent(poly))]
    return evaluate(poly, point)
end


    """
	    func SetParameterValues(ode::ODE{P}, param_values::OrderedDict{P,T}) where {T <: FieldElem,P <: MPolyElem{T}}
    
Substitute parameters with numerical values.

## Input:
- ode, an ODE as above
- param_values, values for (some of) the parameters as dictionary parameter => value

## Output: 
- new ode with the parameters in param_values plugged with the given numbers
    """
function SetParameterValues(ode::ODE{P}, param_values::OrderedDict{P,T}) where {T <: FieldElem,P <: MPolyElem{T}}
    new_vars = map(var_to_str, [v for v in gens(ode.poly_ring) if !(v in keys(param_values))])
    small_ring, small_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_vars)
    eval_dict = OrderedDict(str_to_var(v, ode.poly_ring) => str_to_var(v, small_ring) for v in new_vars)
    merge!(eval_dict, OrderedDict(p => small_ring(val) for (p, val) in param_values))

    return ODE{P}(
        OrderedDict{P,Union{P,Generic.Frac{P}}}(eval_at_dict(v, eval_dict) => eval_at_dict(f, eval_dict) for (v, f) in ode.x_equations),
        OrderedDict{P,Union{P,Generic.Frac{P}}}(eval_at_dict(v, eval_dict) => eval_at_dict(f, eval_dict) for (v, f) in ode.y_equations),
        [eval_at_dict(u, eval_dict) for u in ode.u_vars]
    )
end


    """
        func switch_ring(v::MPolyElem, ring::MPolyRing)

    For a variable v, returns a variable in ring with the same name
    """
function switch_ring(v::MPolyElem, ring::MPolyRing)
    ind = findfirst(vv -> vv == v, gens(parent(v)))
    return str_to_var(string(symbols(parent(v))[ind]), ring)
end


    """
        func var_to_str(v::MPolyElem)

    Convert a variable to type `string`.
    """
function var_to_str(v::MPolyElem)
    ind = findfirst(vv -> vv == v, gens(parent(v)))
    return string(symbols(parent(v))[ind])
end

    """
        func str_to_var(s, ring::MPolyRing)

    Convert a `string`-typed variable to a symbol.
    """
function str_to_var(s, ring::MPolyRing)
    ind = findfirst(v -> (string(v) == s), symbols(ring))
    if ind === nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end


    """
        func unpack_fraction(f::MPolyElem)

    A helper-function, returns a `Tuple` of the input `f` and its parent's multiplicative identity.
    """
function unpack_fraction(f::MPolyElem)
    return (f, one(parent(f)))
end

    """
        func unpack_fraction(f::Generic.Frac{<: MPolyElem})

    A helper-function, returns a `Tuple` of the numerator and denominator of `f`.
    """
function unpack_fraction(f::Generic.Frac{<: MPolyElem})
    return (numerator(f), denominator(f))
end


    """
        func print_for_SIAN(ode::ODE{P}, outputs::Array{P,1}) where P <: MPolyElem{<: FieldElem}

    Prints the ODE in the format accepted by SIAN (https://github.com/pogudingleb/SIAN)
    """
function print_for_SIAN(ode::ODE{P}, outputs::Array{P,1}) where P <: MPolyElem{<: FieldElem}

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

    """
        func macrohelper_extract_vars(equations::Array{Expr,1})

    A helper-function for a macro used in extracting variables from equations.
    """
function macrohelper_extract_vars(equations::Array{Expr,1})
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
# ------------------------------------------------------------------------------

    """
        func macrohelper_clean(ex::Expr)

    A cleanup helper for the macro.
    """
function macrohelper_clean(ex::Expr)
    ex = MacroTools.postwalk(x -> @capture(x, f_'(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> @capture(x, f_(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> x == :(/) ? :(//) : x, ex)
    return ex
end
# ------------------------------------------------------------------------------

    """
    macro ODEmodel(ex::Expr...)

    Macro for creating an ODE from a list of equations and injecting all variables into the global scope.

    Example:

    ```
    ode = @ODEmodel(
        x1'(t) = - a * x1(t),
        y1(t) = x1(t),
    )
    ```
    """
macro ODEmodel(ex::Expr...)
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
    x_dict_create_expr = :($x_dict = OrderedDict{fmpq_mpoly,Union{fmpq_mpoly,Generic.Frac{fmpq_mpoly}}}())
    y_dict_create_expr = :($y_dict = OrderedDict{fmpq_mpoly,Union{fmpq_mpoly,Generic.Frac{fmpq_mpoly}}}())
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


    """
        func generate_replica(ode::ODE{P}, r::Int) where P <: MPolyElem

    Generate a replica of the original input system as per <Theorem here>.
    Returns `ode_r`, and r-fold replica of the original ode.
    States, outputs, and inputs are replicated, parameters are not.
    """
function generate_replica(ode::ODE{P}, r::Int) where P <: MPolyElem
    new_varnames = Array{String}(undef, 0)
    # new_varnames = map(string, ode.parameters)
    for v in vcat(ode.x_vars, ode.y_vars, ode.u_vars)
        append!(new_varnames, [var_to_str(v) * "_r$i" for i in 1:r])
    end
    append!(new_varnames, map(string, ode.parameters))
    new_ring, new_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_varnames)
    new_x_eqs = OrderedDict{P,Union{P,Generic.Frac{P}}}()
    new_y_eqs = OrderedDict{P,Union{P,Generic.Frac{P}}}()
    new_us = Array{P,1}()
    for i in 1:r
        eval = merge(
            OrderedDict(v => str_to_var(var_to_str(v) * "_r$i", new_ring) for v in vcat(ode.x_vars, ode.y_vars, ode.u_vars)),
            OrderedDict(p => switch_ring(p, new_ring) for p in ode.parameters)
        )
        eval_vec = [eval[v] for v in gens(ode.poly_ring)]
        new_x_eqs = merge(
            new_x_eqs,
            OrderedDict{P,Union{P,Generic.Frac{P}}}(evaluate(x, eval_vec) => evaluate(f, eval_vec) for (x, f) in ode.x_equations)
        )
        new_y_eqs = merge(
            new_y_eqs,
            OrderedDict{P,Union{P,Generic.Frac{P}}}(evaluate(x, eval_vec) => evaluate(f, eval_vec) for (x, f) in ode.y_equations)
        )
        append!(new_us, [str_to_var(var_to_str(u) * "_r$i", new_ring) for u in ode.u_vars])
    end
    return ODE{P}(new_x_eqs, new_y_eqs, new_us)
end

##################


    """
        func _reduce_poly_mod_p(poly::MPolyElem{Nemo.fmpq}, p::Int)

    Reduces a polynomial over Q modulo p.
    """
function _reduce_poly_mod_p(poly::MPolyElem{Nemo.fmpq}, p::Int)
    den = denominator(poly)
    num = change_base_ring(Nemo.ZZ, den * poly)
    if Nemo.GF(p)(den) == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
    end
    return change_base_ring(Nemo.GF(p), num) * (1 // Nemo.GF(p)(den))
end

##################


    """
        func make_derivative(var_name, der_order)

    Given a variable name add a derivative order.
    """
function make_derivative(var_name, der_order)
    return(string(var_name, "_", der_order))
end

##################

    """
        func add_to_var(vr, ring, r)

    Convert a variable `vr` to a derivative of order `r` and convert the result to symbol.
    """
function add_to_var(vr, ring, r)
    return str_to_var(make_derivative(vr, r), ring)
end

##################

    """
        func create_jet_ring(var_list, param_list, max_ord)

    Given a list of variables `var_list` and a list of parameters `param_list`, create a jet ring of derivatives up to order `max_ord`.
    """
function create_jet_ring(var_list, param_list, max_ord)
    varnames = vcat(vec(["$(s)_$i" for s in var_list, i in 0:max_ord]), "z_aux", ["$(s)_0" for s in param_list])
    return Nemo.PolynomialRing(Nemo.QQ, varnames)[1]
end

##################

    """
        func differentiate_all(diff_poly, var_list, shft, max_ord)

    Differentiate a polynomial `diff_poly` with respect to `var_list` up to `max_ord` order.
    """
function differentiate_all(diff_poly, var_list, shft, max_ord)
    result = 0
    for i in 1:(shft * (max_ord + 1))
        result = result + derivative(diff_poly, var_list[i]) * var_list[i + shft]
    end
    return result
end

##################

# function unpack_fraction(f)
#    if applicable(numerator, f)
#        return (numerator(f), denominator(f))
#    end
#    return (f, parent(f)(1))
# end

##################

    """
        func eval_frac(frac, vars, vals)

    Evaluate a given fraction `frac` with values `vals` in place of variables `vars`.
    """
function eval_frac(frac, vars, vals)
    fr = unpack_fraction(frac)
    return (evaluate(fr[1], vars, vals) // evaluate(fr[2], vars, vals))
end

###################


    """
        func sample_point(bound, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q)

    Sample random values for parameters of the polynomial system.
    """
function sample_point(bound, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q)
    local u_hat, theta_hat, all_subs

    s = length(all_params)
    y_hat_vars = Array{fmpq_mpoly}(undef, 0)
    y_hat_vals = Array{fmpq}(undef, 0)

    while true
        theta_hat = [fmpq(rnd) for rnd in rand(0:bound, s)]
        u_hat =  [fmpq(rnd) for rnd in rand(0:bound, length(u_variables))]
        all_subs = [vcat(all_params, u_variables), vcat(theta_hat, u_hat)]
        if evaluate(Q, all_subs[1], all_subs[2]) == 0
            next
        end
        vls = insert_zeros_to_vals(all_subs[1], all_subs[2])
        for i in 0:(s + 1)
            for j in 1:length(y_vars)
                eq = Y_eq[(j - 1) * (s + 2) + i + 1][2]
                vl = evaluate(unpack_fraction(eq)[1], vls) // evaluate(unpack_fraction(eq)[2], vls)
                y_hat_vars = vcat(y_hat_vars, Y_eq[(j - 1) * (s + 2) + i + 1][1])
                y_hat_vals = vcat(y_hat_vals, vl)
                vls[var_index(Y_eq[(j - 1) * (s + 2) + i + 1][1])] = vl
                all_subs = [vcat(all_subs[1], Y_eq[(j - 1) * (s + 2) + i + 1][1]), vcat(all_subs[2], vl)]
            end
            for j in 1:length(x_vars)
                eq = X_eq[(j - 1) * (s + 2) + i + 1][2]
                vl = evaluate(unpack_fraction(eq)[1], vls) // evaluate(unpack_fraction(eq)[2], vls)
                all_subs = [vcat(all_subs[1], X_eq[(j - 1) * (s + 2) + i + 1][1]),vcat(all_subs[2], vl)]
                vls[var_index(X_eq[(j - 1) * (s + 2) + i + 1][1])] = vl
            end
        end
        break
    end
    return([[y_hat_vars, y_hat_vals], [u_variables, u_hat], [all_params, theta_hat], all_subs])
end

##################
    """
    func jacobi_matrix(pol_arr, vrs, vals)

    Generate a Jacobi matrix from a given array of polynomial `pol_arr`,
    with respect to variables `vars`.
    The matrix is evaluated at `vals` from a symbolic to numeric representation.
    """
function jacobi_matrix(pol_arr, vrs, vals)
    m = Nemo.MatrixSpace(Nemo.QQ, length(pol_arr), length(vrs))()
    for i in 1:length(pol_arr)
        for j in 1:length(vrs)
            m[i, j] = evaluate(derivative(pol_arr[i], vrs[j]), vals)
        end
    end
    return m
end

    ##################
    """
        func get_order_var2(diff_var, non_jet_vars, shft, s)
    """
function get_order_var2(diff_var, non_jet_vars, shft, s)
    idx = var_index(diff_var)
    if idx <= shft * (s + 3)
        return([non_jet_vars[rem(idx - 1, shft) + 1], div(idx - 1, shft)])
    else
        return([non_jet_vars[idx - shft * (s + 2) - 1], 0])
    end
end

##################

    """
        func get_order_var(diff_var, non_jet_ring)

        A helper function to obtain derivative order from string.
    """
function get_order_var(diff_var, non_jet_ring)
    rex = match(r"^(.*_)([0-9]+)$", string(diff_var))
    if rex === nothing
        return(["", ""])
    else
        return([str_to_var(first(rex[1], length(rex[1]) - 1), non_jet_ring), parse(Int, rex[2])])
    end
end

##################

    """
        func get_vars(diff_poly, var_list, non_jet_vars, shft, s)

        Get variables from `diff_poly` based on the intersection with `var_list`.
    """
function get_vars(diff_poly, var_list, non_jet_vars, shft, s)
    return [v for v in vars(diff_poly) if get_order_var2(v, non_jet_vars, shft, s)[1] in var_list]
end


##################


    """
        func compare_diff_var(dvl, dvr, non_jet_vars, shft, s)
    
        Comparison method of variables based on order.
    """
function compare_diff_var(dvl, dvr, non_jet_vars, shft, s)
    vl, hl = get_order_var2(dvl, non_jet_vars, shft, s)
    vr, hr = get_order_var2(dvr, non_jet_vars, shft, s)
    if hl != hr
        return (hl > hr)
    end
    if length(string(vl)) != length(string(vr))
        return (length(string(vl)) > length(string(vr)))
    end
    return (vr >= vl)
end

    ##################
    """
        func parent_ring_change(poly::MPolyElem, new_ring::MPolyRing)

    Converts a polynomial to a different polynomial ring.

    ## Input:
      - `poly::MPolyElem` - a polynomial to be converted
      - `new_ring::MPolyRing` - a polynomial ring such that every variable name
          appearing in poly appears among the generators
    
    ## Output:
      - a polynomial in new_ring "equal" to `poly`
    """
function parent_ring_change(poly::MPolyElem, new_ring::MPolyRing)

    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = Array{Any,1}()
        
    for u in symbols(old_ring)
        push!(
            var_mapping,
            findfirst(v -> (string(u) == string(v)), symbols(new_ring))
        )
    end
    builder = MPolyBuildCtx(new_ring)
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        new_exp = [0 for _ in gens(new_ring)]
        for i in 1:length(exp)
            if exp[i] != 0
                if var_mapping[i] === nothing
                    throw(Base.ArgumentError("The polynomial contains a variable not present in the new ring $poly"))
                else
                    new_exp[var_mapping[i]] = exp[i]
        end
            end
        end
        push_term!(builder, new_ring.base_ring(coef), new_exp)
    end
    return finish(builder)
end

    ##################
    """
        func insert_zeros_to_vals(var_arr, val_arr)

    Insert zeros at positions based on the variables' index.

    """
function insert_zeros_to_vals(var_arr, val_arr)
    ### worstcase O(n^2) ??
    all_val_arr = zeros(fmpq, length(gens(parent(var_arr[1]))))
    for i in 1:length(var_arr)
        all_val_arr[var_index(var_arr[i])] = val_arr[i]
    end
    return all_val_arr
end

    ##################
    """
        func add_zero_to_vars(poly::MPolyElem, new_ring::MPolyRing)
Converts a polynomial to a different polynomial ring.

## Input

  - `poly::MPolyElem` - a polynomial to be converted
  - `new_ring::MPolyRing` - a polynomial ring such that every variable name
      appearing in poly appears among the generators

## Output
    -  a polynomial in new_ring "equal" to `poly`
    """
function add_zero_to_vars(poly::MPolyElem, new_ring::MPolyRing)
    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = Array{Any,1}()
    for u in symbols(old_ring)
        push!(
            var_mapping,
            findfirst(v -> (string(u, "_0") == string(v)), symbols(new_ring))
        )
    end
    builder = MPolyBuildCtx(new_ring)
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        new_exp = [0 for _ in gens(new_ring)]
        for i in 1:length(exp)
            if exp[i] != 0
                if var_mapping[i] === nothing
                    throw(Base.ArgumentError("The polynomial contains a variable not present in the new ring $poly"))
                else
                    new_exp[var_mapping[i]] = exp[i]
        end
            end
        end
        push_term!(builder, new_ring.base_ring(coef), new_exp)
    end
    return finish(builder)
end

#############

    """ 
    function get_parameters(ode; initial_conditions=true)
Retrieve parameters from the `ode` system. Retrieve initial conditions if `initial_conditions` is set `true`.

## Input
    - `ode::ODE` - an ODE system
    - `initial_conditions::Bool` - whether to extract initial conditions. Default `true`.

## Output
        - Array of parameters (and initial conditions).
    """
function get_parameters(ode; initial_conditions=true)
    if initial_conditions
        return vcat(ode.parameters, ode.x_vars)
    else
        return ode.parameters
    end
end
        
    """
        func var_to_symb(var)

    Convert a variable `var` to `symbol`.
    """
function var_to_symb(var)
    symbols(parent(var))[var_index(var)]
end

    """
        func add_to_vars_in_replica(poly::MPolyElem, mu, new_ring::MPolyRing, r)

    A helper routine to add variables from symbols of the old ring based on `poly`, 
    to the `new_ring` object.
    """
function add_to_vars_in_replica(poly::MPolyElem, mu, new_ring::MPolyRing, r)
    
    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = Array{Any,1}()
    mu_symbols = [var_to_symb(m) for m in mu]
        
    for u in symbols(old_ring)
        if u in mu_symbols
            push!(
                var_mapping,
                findfirst(v -> (string(u) == string(v)), symbols(new_ring))
            )
        else
            push!(
                var_mapping,
                findfirst(v -> (string(u, "_", r) == string(v)), symbols(new_ring))
            )
        end
    end
    builder = MPolyBuildCtx(new_ring)
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        new_exp = [0 for _ in gens(new_ring)]
        for i in 1:length(exp)
            if exp[i] != 0
                if var_mapping[i] === nothing
                    throw(Base.ArgumentError("The polynomial contains a variable not present in the new ring $poly"))
                else
                    new_exp[var_mapping[i]] = exp[i]
                end
            end
        end
        push_term!(builder, new_ring.base_ring(coef), new_exp)
    end
    return finish(builder)
end


    ############# Main Code
    """
        func identifiability_ode(ode, params_to_assess; p=0.99, p_mod=0, nthrds=64)
    
    Perform identifiability check for a given `ode` system with respect to parameters in `params_to_assess` list.
    
## Input

  - `ode` - an ODE system returned by the `@ODEmodel` macro.
  - `params_to_assess` - an array of parameters returned by `get_parameters` function.
  - `p` - probability of correctness, default `0.99`.
  - `p_mod` - a prime characteristic, default `0`.
    - `nthrds` - number of threads for concurrency, default `64`.
    """
function identifiability_ode(ode, params_to_assess; p=0.99, p_mod=0, nthrds=64)

    println("Solving the problem")
    # 1.Construct the maximal system
    
    # (a) ---------------
    println("Constructing the maximal system")

    non_jet_ring = ode.poly_ring
    all_indets = gens(non_jet_ring)
    x_vars = ode.x_vars
    y_vars = ode.y_vars
    u_vars = ode.u_vars
    mu = ode.parameters
    p_local = p + length(params_to_assess) * 10^(-18)

    n = length(x_vars)
    m = length(y_vars)
    u = length(u_vars)
    s = length(mu) + n

    Rjet = create_jet_ring(vcat(x_vars, y_vars, u_vars), mu, s + 2)
    gens_Rjet = gens(Rjet)
    z_aux = gens_Rjet[end - length(mu)]
    x_eqs = collect(values(ode.x_equations))
    y_eqs = collect(values(ode.y_equations))

    x_eqs = [[add_to_var(x_vars[i], Rjet, 1), add_zero_to_vars(unpack_fraction(x_eqs[i])[1], Rjet) // add_zero_to_vars(unpack_fraction(x_eqs[i])[2], Rjet)] for i in 1:n]
    y_eqs = [[add_to_var(y_vars[i], Rjet, 0), add_zero_to_vars(unpack_fraction(y_eqs[i])[1], Rjet) // add_zero_to_vars(unpack_fraction(y_eqs[i])[2], Rjet)] for i in 1:m]

    eqs = vcat(x_eqs, y_eqs)
    Q = foldl(lcm, [unpack_fraction(ex[2])[2] for ex in eqs])

    not_int_cond_params = gens_Rjet[(end - length(ode.parameters) + 1):end]
    all_params = vcat(not_int_cond_params, gens_Rjet[1:n])
    x_variables = gens_Rjet[1:n]
    for i in 1:(s + 1)
        x_variables = vcat(x_variables, gens_Rjet[(i * (n + m + u) + 1):(i * (n + m + u) + n)])
    end
    u_variables = gens_Rjet[(n + m + 1):(n + m + u)]
    for i in 1:(s + 1)
        u_variables = vcat(u_variables, gens_Rjet[((n + m + u) * i + n + m + 1):((n + m + u) * (i + 1))])
    end

    # (b,c) -------------
    X = Array{fmpq_poly}(undef, 0)
    X_eq = Array{fmpq_poly}(undef, 0)
    for i in 1:n
        X = vcat(X, [Array{fmpq_poly}(undef, 0)])
        poly_d = unpack_fraction(x_eqs[i][1] - x_eqs[i][2])[1]
        for j in 0:s + 1
            if j > 0
                poly_d = differentiate_all(poly_d, gens_Rjet, n + m + u, j)
            end
            leader = gens_Rjet[i + (n + m + u) * (j + 1)]
            separant = derivative(poly_d, leader)
            X[i] = vcat(X[i], poly_d)
            X_eq = vcat(X_eq, [[leader,-(poly_d - separant * leader) // separant]])
        end
    end

    # (d,e) -----------
    Y = Array{fmpq_poly}(undef, 0)
    Y_eq = Array{fmpq_poly}(undef, 0)
    for i in 1:m
        Y = vcat(Y, [Array{fmpq_poly}(undef, 0)])
        poly_d = unpack_fraction(y_eqs[i][1] - y_eqs[i][2])[1]
        for j in 0:s + 1
            if j > 0
                poly_d = differentiate_all(poly_d, gens_Rjet, n + m + u, j - 1)
            end
            leader = gens_Rjet[i + n + (n + m + u) * j]
            separant = derivative(poly_d, leader)
            Y[i] = vcat(Y[i], poly_d)
            Y_eq = vcat(Y_eq, [[leader,-(poly_d - separant * leader) // separant]])
        end
    end

    # 2.Truncate

    println("Truncating")

    # (a) -----------------------
    
    d0 = BigInt(maximum(vcat([total_degree(unpack_fraction(Q * eq[2])[1]) for eq in eqs], total_degree(Q))))
    
    # (b) -----------------------
    
    D1 = floor(BigInt, (length(params_to_assess) + 1) * 2 * d0 * s * (n + 1) * (1 + 2 * d0 * s) / (1 - p_local))

    # (c, d) --------------------
    
    sample = sample_point(D1, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q)
    all_subs = sample[4]
    u_hat = sample[2]
    y_hat = sample[1]

    # (e) -----------------------

    alpha = [1 for i in 1:n]
    beta = [0 for i in 1:m]
    Et = Array{fmpq_poly}(undef, 0)
    x_theta_vars = all_params
    prolongation_possible = [1 for i in 1:m]
    
    # (f) -----------------------
    all_x_theta_vars_subs = insert_zeros_to_vals(all_subs[1], all_subs[2])
    eqs_i_old = Array{fmpq_mpoly}(undef, 0)
    evl_old = Array{fmpq_mpoly}(undef, 0)
    while sum(prolongation_possible) > 0
        for i in 1:m
            if prolongation_possible[i] == 1
                eqs_i = vcat(Et, Y[i][beta[i] + 1])
                evl = [evaluate(eq, vcat(u_hat[1], y_hat[1]), vcat(u_hat[2], y_hat[2])) for eq in eqs_i if !(eq in eqs_i_old)]
                evl_old = vcat(evl_old, evl)
                JacX = jacobi_matrix(evl_old, x_theta_vars, all_x_theta_vars_subs)
                eqs_i_old = eqs_i
                if LinearAlgebra.rank(JacX) == length(eqs_i)
                    Et = vcat(Et, Y[i][beta[i] + 1])
                    beta[i] = beta[i] + 1
                    # adding necessary X-equations
                    polys_to_process = vcat(Et, [Y[k][beta[k] + 1] for k in 1:m])
                    while length(polys_to_process) != 0
                        new_to_process = Array{fmpq_mpoly}(undef, 0)
                        vrs = Set{fmpq_mpoly}()
                        for poly in polys_to_process
                            vrs = union(vrs, [v for v in vars(poly) if v in x_variables])
                        end
                        vars_to_add = Set{fmpq_mpoly}(v for v in vrs if !(v in x_theta_vars))
                        for v in vars_to_add
                            x_theta_vars = vcat(x_theta_vars, v)
                            ord_var = get_order_var2(v, all_indets, n + m + u, s)
                            var_idx = var_index(ord_var[1])
                            poly = X[ var_idx ][ ord_var[2] ]
                            Et = vcat(Et, poly)
                            new_to_process = vcat(new_to_process, poly)
                            alpha[ var_idx ] = max(alpha[ var_idx ], ord_var[2] + 1)
                        end
                        polys_to_process = new_to_process
                    end
                else
                    prolongation_possible[i] = 0
                end
            end
        end
    end

    println("Assessing local identifiability")

    max_rank = length(Et)
    for i in 1:m
        for j in (beta[i] + 1):length(Y[i])
            to_add = true
            for v in get_vars(Y[i][j], x_vars, all_indets, n + m + u, s)
                if !(v in x_theta_vars)
                    to_add = false
                end
            end
            if to_add
                beta[i] = beta[i] + 1
                Et = vcat(Et, Y[i][j])
            end
        end
    end

    theta_l = Array{fmpq_mpoly}(undef, 0)
    params_to_assess = [add_to_var(param, Rjet, 0) for param in params_to_assess]
    Et_eval_base = [evaluate(e, vcat(u_hat[1], y_hat[1]), vcat(u_hat[2], y_hat[2])) for e in Et]
        for param_0 in params_to_assess
        other_params = [v for v in x_theta_vars if v != param_0]
        Et_subs = [evaluate(e, [param_0], [evaluate(param_0, all_x_theta_vars_subs)]) for e in Et_eval_base]
        JacX = jacobi_matrix(Et_subs, other_params, all_x_theta_vars_subs)
        if LinearAlgebra.rank(JacX) != max_rank
            theta_l = vcat(theta_l, param_0)
        end
    end
    if length(theta_l) == 0
        println("\n=== Summary ===")
        println("Globally identifiable parameters:                 []")
        println("Locally but not globally identifiable parameters: []")
        println("Not identifiable parameters:                      [", join(params_to_assess, ", "), "]")
    else
        println("Locally identifiable parameters: [", join([get_order_var(th, non_jet_ring)[1] for th in theta_l], ", "), "]")
        println("Not identifiable parameters:     [", join([get_order_var(th, non_jet_ring)[1] for th in setdiff(params_to_assess, theta_l)], ", "), "]")
    end
    # 3. Randomize.

    println("Randomizing")
    # (a) ------------
    deg_variety =  foldl(*, [BigInt(total_degree(e)) for e in Et])
    D2 = floor(BigInt, 6 * length(theta_l) * deg_variety * (1 + 2 * d0 * maximum(beta)) / (1 - p_local))
    # (b, c) ---------
    sample = sample_point(D2, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q)
    y_hat = sample[1]
    u_hat = sample[2]
    theta_hat = sample[3]
    
    # (d) ------------
    Et_hat = [evaluate(e, vcat(y_hat[1], u_hat[1]), vcat(y_hat[2], u_hat[2])) for e in Et]
    Et_x_vars = Set{fmpq_mpoly}()
    for poly in Et_hat
        Et_x_vars = union(Et_x_vars, get_vars(poly, x_vars, all_indets, n + m + u, s))
    end
    Q_hat = evaluate(Q, u_hat[1], u_hat[2])
    vrs_sorted = vcat(sort([e for e in Et_x_vars], lt=(x, y) -> compare_diff_var(x, y, all_indets, n + m + u, s)), z_aux, sort(not_int_cond_params, rev=true))
    # 4. Determine.
    println("GB computation")


    if p_mod > 0
        Et_hat = [_reduce_poly_mod_p(e, p_mod) for e in Et_hat]
        z_aux = _reduce_poly_mod_p(z_aux, p_mod)
        Q_hat = _reduce_poly_mod_p(Q_hat, p_mod)
        Rjet_new, vrs_sorted = Singular.PolynomialRing(Singular.Fp(p_mod), [string(v) for v in vrs_sorted], ordering=:degrevlex)
    else
        Rjet_new, vrs_sorted = Singular.PolynomialRing(Singular.QQ, [string(v) for v in vrs_sorted], ordering=:degrevlex)
    end

    theta_g = Array{spoly}(undef, 0)
    Et_hat = [parent_ring_change(e, Rjet_new) for e in Et_hat]
    gb = GroebnerBasis.f4(Ideal(Rjet_new, vcat(Et_hat, parent_ring_change(z_aux * Q_hat, Rjet_new) - 1)), nthrds=nthrds)
    println("Remainder computation")

    if p_mod > 0
        theta_l_new = [parent_ring_change(_reduce_poly_mod_p(th, p_mod), Rjet_new) for th in theta_l]
            for i in 1:length(theta_l)
            if Singular.reduce(theta_l_new[i], gb) == parent_ring_change(_reduce_poly_mod_p(Rjet(theta_hat[2][findfirst(isequal(theta_l[i]), theta_hat[1])]), p_mod), Rjet_new)
                theta_g = vcat(theta_g, theta_l_new[i])
            end
        end
    else
        theta_l_new = [parent_ring_change(th, Rjet_new) for th in theta_l]
            for i in 1:length(theta_l)
            if Singular.reduce(theta_l_new[i], gb) == parent_ring_change(Rjet(theta_hat[2][findfirst(isequal(theta_l[i]), theta_hat[1])]), Rjet_new)
                theta_g = vcat(theta_g, theta_l_new[i])
            end
        end
    end
    println("\n=== Summary ===")
    println("Globally identifiable parameters:                 [", join([get_order_var(th, non_jet_ring)[1] for th in theta_g], ", "), "]")
    println("Locally but not globally identifiable parameters: [", join([get_order_var(th, non_jet_ring)[1] for th in setdiff(theta_l_new, theta_g)], ", "), "]")
    println("Not identifiable parameters:                      [", join([get_order_var(th, non_jet_ring)[1] for th in setdiff(params_to_assess, theta_l)], ", "), "]")
    println("===============")
end

end # module SIAN
