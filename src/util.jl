# ------------------------------------------------------------------------------
"""
    func eval_at_dict(poly::P, d::OrderedDict{P,<: RingElem}) where P <: MPolyElem

Evaluates a polynomial on a dict var => val
missing values are replaced with zeroes
"""
function eval_at_dict(poly::P, d::OrderedDict{P,<: RingElem}) where P <: MPolyElem
    
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
function unpack_fraction(f::Generic.Frac{<: MPolyElem})
    return (numerator(f), denominator(f))
end

# ------------------------------------------------------------------------------
"""
    func eval_frac(frac, vars, vals)

Evaluate a given fraction `frac` with values `vals` in place of variables `vars`.
"""
function eval_frac(frac, vars, vals)
    fr = unpack_fraction(frac)
    return(evaluate(fr[1], vars, vals) // evaluate(fr[2], vars, vals))
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
        push_term!(builder, new_ring.base_ring(coef), new_exp)
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
