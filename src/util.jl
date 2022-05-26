# Comment: all the function except for
#   - insert_zeros_to_vals
#   - add_zero_to_vars
#   - var_to_symb
#   - add_to_vars_in_replica
#   Have been adapted from a different repository. A precise reference will be
#   inserted once that repository will become public


# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
""" 
    func get_parameters(ode; initial_conditions=true)

Retrieve parameters from the `ode` system. Retrieve initial conditions if `initial_conditions` is set `true`.

## Input
    - `ode::ODE` - an ODE system
    - `initial_conditions::Bool` - whether to extract initial conditions. Default `true`.

## Output
        - Array of parameters (and initial conditions).
"""
function get_parameters(ode; initial_conditions = true)
    if initial_conditions
        return vcat(ode.parameters, ode.x_vars)
    else
        return ode.parameters
    end
end


# ------------------------------------------------------------------------------
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
    return (m)
end

# ------------------------------------------------------------------------------
"""
    func get_order_var2(diff_var, non_jet_vars, shft, s)
"""
function get_order_var2(diff_var, non_jet_vars, shft, s)
    idx = var_index(diff_var)
    if idx <= shft * (s + 3)
        return ([non_jet_vars[rem(idx - 1, shft)+1], div(idx - 1, shft)])
    else
        return ([non_jet_vars[idx-shft*(s+2)-1], 0])
    end
end

# ------------------------------------------------------------------------------
"""
    func get_order_var(diff_var, non_jet_ring)

A helper function to obtain derivative order from string.
"""
function get_order_var(diff_var, non_jet_ring)
    rex = match(r"^(.*_)([0-9]+)$", string(diff_var))
    if rex === nothing
        return (["", ""])
    else
        return ([str_to_var(first(rex[1], length(rex[1]) - 1), non_jet_ring), parse(Int, rex[2])])
    end
end

# ------------------------------------------------------------------------------
"""
    func get_vars(diff_poly, var_list, non_jet_vars, shft, s)

Get variables from `diff_poly` based on the intersection with `var_list`.
"""
function get_vars(diff_poly, var_list, non_jet_vars, shft, s)
    return [v for v in vars(diff_poly) if get_order_var2(v, non_jet_vars, shft, s)[1] in var_list]
end


"""
    func make_derivative(var_name, der_order)

Given a variable name `var_name` add a derivative order `der_order`.
"""
function make_derivative(var_name, der_order)
    return (string(var_name, "_", der_order))
end

# ------------------------------------------------------------------------------
"""
    func add_to_var(vr, ring, r)

Convert a variable `vr` to a derivative of order `r` and convert the result to symbol.
"""
function add_to_var(vr, ring, r)
    return str_to_var(make_derivative(vr, r), ring)
end

# ------------------------------------------------------------------------------
"""
    func create_jet_ring(var_list, param_list, max_ord)

Given a list of variables `var_list` and a list of parameters `param_list`, create a jet ring of derivatives up to order `max_ord`.
"""
function create_jet_ring(var_list, param_list, max_ord)
    varnames = vcat(vec(["$(s)_$i" for s in var_list, i in 0:max_ord]), "z_aux", ["$(s)_0" for s in param_list])
    return Nemo.PolynomialRing(Nemo.QQ, varnames)[1]
end

# ------------------------------------------------------------------------------
"""
    func differentiate_all(diff_poly, var_list, shft, max_ord)

Differentiate a polynomial `diff_poly` with respect to `var_list` up to `max_ord` order.
"""
function differentiate_all(diff_poly, var_list, shft, max_ord)
    result = 0
    for i in 1:(shft*(max_ord+1))
        result = result + derivative(diff_poly, var_list[i]) * var_list[i+shft]
    end
    return (result)
end

# ------------------------------------------------------------------------------
"""
    func eval_at_dict(poly::P, d::OrderedDict{P,<: RingElem}) where P <: MPolyElem

Evaluates a polynomial on a dict var => val
missing values are replaced with zeroes
"""
function eval_at_dict(poly::P, d::OrderedDict{P,<:RingElem}) where {P<:MPolyElem}

    point = [get(d, v, base_ring(parent(poly))(0)) for v in gens(parent(poly))]
    return evaluate(poly, point)
end

# ------------------------------------------------------------------------------
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
    if ind == nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end

# ------------------------------------------------------------------------------

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
function unpack_fraction(f::Generic.Frac{<:MPolyElem})
    return (numerator(f), denominator(f))
end

# ------------------------------------------------------------------------------
"""
    func eval_frac(frac, vars, vals)

Evaluate a given fraction `frac` with values `vals` in place of variables `vars`.
"""
function eval_frac(frac, vars, vals)
    fr = unpack_fraction(frac)
    return (evaluate(fr[1], vars, vals) // evaluate(fr[2], vars, vals))
end

# ------------------------------------------------------------------------------
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
                if var_mapping[i] == nothing
                    throw(Base.ArgumentError("The polynomial contains a variable not present in the new ring $poly"))
                else
                    new_exp[var_mapping[i]] = exp[i]
                end
            end
        end
        if typeof(coef) <: Nemo.fmpq
            push_term!(builder, new_ring.base_ring(coef), new_exp)
        else
            push_term!(builder, new_ring.base_ring(Nemo.data(coef)), new_exp)
        end
    end
    return finish(builder)
end

# ------------------------------------------------------------------------------
"""
    func insert_zeros_to_vals(var_arr, val_arr)

Insert zeros at positions based on the variables' index.
"""
function insert_zeros_to_vals(var_arr, val_arr)
    all_val_arr = zeros(fmpq, length(gens(parent(var_arr[1]))))
    for i in 1:length(var_arr)
        all_val_arr[var_index(var_arr[i])] = val_arr[i]
    end
    return all_val_arr
end

# ------------------------------------------------------------------------------
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
                if var_mapping[i] == nothing
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

# ------------------------------------------------------------------------------
"""
    func var_to_symb(var)

Convert a variable `var` to `symbol`.
"""
function var_to_symb(gn)
    symbols(parent(gn))[var_index(gn)]
end

# ------------------------------------------------------------------------------
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
                if var_mapping[i] == nothing
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

# ------------------------------------------------------------------------------

