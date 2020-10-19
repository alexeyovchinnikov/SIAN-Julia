println("Loading packages")
using Oscar
using LinearAlgebra
using Singular
using GroebnerBasis
using Dates

##################

function _reduce_poly_mod_p(poly::MPolyElem{Nemo.fmpq}, p::Int)
    """
    Reduces a polynomial over Q modulo p
    """
    den = denominator(poly)
    num = change_base_ring(Nemo.ZZ, den * poly)
    if GF(p)(den) == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
    end
    return change_base_ring(Fp(p), num) * (1 // Fp(p)(den))
end

##################

function str_to_var2(s, ring)
    gns = gens(ring)
    return gns[findfirst(v -> (string(v) == s), gns)]
end


function str_to_var(s, ring::MPolyRing)
    ind = findfirst(v -> (string(v) == s), symbols(ring))
    if ind == nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end

##################

function make_derivative(var_name, der_order)
    return(string(var_name, "_", der_order))
end

##################

function add_to_var(vr, ring, r)
    return str_to_var(make_derivative(vr, r), ring)
end

##################

function create_jet_ring(var_list, param_list, max_ord)
    varnames = vcat(vec(["$(s)_$i" for s in var_list, i in 0:max_ord]), "z_aux", ["$(s)_0" for s in param_list])
    return Nemo.PolynomialRing(Nemo.QQ, varnames)[1]
end

##################

function differentiate_all(diff_poly, var_list, shft, max_ord)
    result = 0
    for i in 1:(shft * (max_ord + 1))
        result = result + derivative(diff_poly, var_list[i]) * var_list[i+shft]
    end
    return(result)
end

##################

function unpack_fraction(f)
    if applicable(numerator, f)
        return (numerator(f), denominator(f))
    end
    return (f, parent(f)(1))
end

##################

function eval_frac(frac, vars, vals)
    fr = unpack_fraction(frac)
    return(evaluate(fr[1], vars, vals)//evaluate(fr[2], vars, vals))
end

###################

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
                vls[var_index(Y_eq[(j - 1)*(s + 2) + i + 1][1])] = vl
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

function jacobi_matrix(pol_arr, vrs, vals)
    m = Nemo.MatrixSpace(Nemo.QQ, length(pol_arr), length(vrs))()
    for i in 1:length(pol_arr)
        for j in 1:length(vrs)
            m[i, j] = evaluate(derivative(pol_arr[i], vrs[j]), vals)
        end
    end
    return(m)
end

##################

function get_order_var2(diff_var, non_jet_vars, shft, s) 
    idx = var_index(diff_var)
    if idx <= shft * (s + 3)
        return([non_jet_vars[rem(idx - 1, shft) + 1], div(idx - 1, shft)])
    else
        return([non_jet_vars[idx - shft*(s + 2) - 1], 0])
    end
end

##################

function get_order_var(diff_var, non_jet_ring)
    rex = match(r"^(.*_)([0-9]+)$", string(diff_var))
    if rex === nothing
        return(["", ""])
    else
        return([str_to_var(first(rex[1], length(rex[1])-1), non_jet_ring), parse(Int, rex[2])])
    end
end

##################

function get_vars(diff_poly, var_list, non_jet_vars, shft, s)
    return [v for v in vars(diff_poly) if get_order_var2(v, non_jet_vars, shft, s)[1] in var_list]
end


##################

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

function parent_ring_change(poly::MPolyElem, new_ring::MPolyRing)
    """
    Converts a polynomial to a different polynomial ring
    Input
      - poly - a polynomial to be converted
      - new_ring - a polynomial ring such that every variable name
          appearing in poly appears among the generators
    Output: a polynomial in new_ring "equal" to poly
    """
    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = Array{Any, 1}()
    
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

##################

function insert_zeros_to_vals(var_arr, val_arr)
    all_val_arr = zeros(fmpq, length(gens(parent(var_arr[1]))))
    for i in 1:length(var_arr)
        all_val_arr[var_index(var_arr[i])] = val_arr[i]
    end 
    return all_val_arr
end

##################

function add_zero_to_vars(poly::MPolyElem, new_ring::MPolyRing)
    """
    Converts a polynomial to a different polynomial ring
    Input
      - poly - a polynomial to be converted
      - new_ring - a polynomial ring such that every variable name
          appearing in poly appears among the generators
    Output: a polynomial in new_ring "equal" to poly
    """
    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = Array{Any, 1}()
    for u in symbols(old_ring)
        push!(
            var_mapping,
            findfirst(v -> (string(u,"_0") == string(v)), symbols(new_ring))
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

#############

function get_parameters(x_eqs, y_eqs, u_vars; initial_conditions = true)
    x_vars = [xv[1] for xv in x_eqs]
    y_vars = [yv[1] for yv in y_eqs]
    all_vars = vcat(x_vars, y_vars, u_vars)
    
    all_indets = vars(x_eqs[1][2])
    for i in 2:length(x_eqs)
        all_indets = union(all_indets, vars(unpack_fraction(x_eqs[i][2])[1]), vars(unpack_fraction(x_eqs[i][2])[2]))
    end

    for i in 1:length(y_eqs)
        all_indets = union(all_indets, vars(unpack_fraction(y_eqs[i][2])[1]), vars(unpack_fraction(y_eqs[i][2])[2]))
    end

    #finding which ones are the parameters
    if initial_conditions
        return vcat(setdiff(all_indets, all_vars), x_vars)
    else
        return setdiff(all_indets, all_vars)
    end 
end

function var_to_symb(gn)
   symbols(parent(gn))[var_index(gn)]
end

function add_to_vars_in_replica(poly::MPolyElem, mu, new_ring::MPolyRing, r)
    """
    Converts a polynomial to a different polynomial ring
    Input
      - poly - a polynomial to be converted
      - new_ring - a polynomial ring such that every variable name
          appearing in poly appears among the generators
    Output: a polynomial in new_ring "equal" to poly
    """
    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = Array{Any, 1}()
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
                findfirst(v -> (string(u,"_",r) == string(v)), symbols(new_ring))
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


function generate_replica(x_eqs, y_eqs, u_vars, r)
    non_jet_ring = parent(x_eqs[1][1])
    x_vars = [xv[1] for xv in x_eqs]
    y_vars = [yv[1] for yv in y_eqs]   
    mu = setdiff(get_parameters(x_eqs, y_eqs, u_vars), x_vars)
    varnames = vcat(vec(["$(s)_$i" for s in vcat(x_vars, y_vars, u_vars), i in 1:r]),  ["$(s)" for s in mu])
    replica_ring = Nemo.PolynomialRing(Nemo.QQ, varnames)[1]
    x_eqs_replica = Array{Array{fmpq_mpoly,fmpq_mpoly}}(undef, 0)
    y_eqs_replica = Array{Array{fmpq_mpoly,fmpq_mpoly}}(undef, 0)
    u_vars_replica = Array{fmpq_mpoly}(undef, 0)
    for i in 1:r 
        x_eqs_replica = vcat(x_eqs_replica, [[add_to_var(x_vars[j], replica_ring, i), add_to_vars_in_replica(x_eqs[j][2], mu, replica_ring, i)] for j in 1:length(x_eqs)]) 
        y_eqs_replica = vcat(y_eqs_replica, [[add_to_var(y_vars[j], replica_ring, i), add_to_vars_in_replica(y_eqs[j][2], mu, replica_ring, i)] for j in 1:length(y_eqs)])
        u_vars_replica = vcat(u_vars_replica, [add_to_var(u_var, replica_ring, i) for u_var in u_vars])
    end
    return(x_eqs_replica, y_eqs_replica, u_vars_replica)
end

############# Main Code

function identifiability_ode(x_eqs, y_eqs, u_vars, params_to_assess, p)

    println("Solving the problem")
# 1.Construct the maximal system
    
    # (a) ---------------
    non_jet_ring_old = parent(x_eqs[1][1])
    x_vars = [xv[1] for xv in x_eqs]
    y_vars = [yv[1] for yv in y_eqs]
    println("Constructing the maximal system")
    
    all_vars   = vcat(x_vars, y_vars, u_vars)
    
    #computing all symbols in the equations
    all_indets = vars(x_eqs[1][2])
    for i in 2:length(x_eqs)
        all_indets = union(all_indets, vars(unpack_fraction(x_eqs[i][2])[1]), vars(unpack_fraction(x_eqs[i][2])[2]))
    end
    for i in 1:length(y_eqs)
        all_indets = union(all_indets, vars(unpack_fraction(y_eqs[i][2])[1]), vars(unpack_fraction(y_eqs[i][2])[2]))
    end
    
    #finding which ones are the parameters
    mu = setdiff(all_indets, all_vars)
    all_indets = vcat(all_vars, mu)
    

    non_jet_ring, all_indets = Nemo.PolynomialRing(Nemo.QQ, ["$(s)" for s in all_indets])
    
    p_local = p + length(params_to_assess) * 10^(-18)

    n = length(x_vars)
    m = length(y_vars)
    u = length(u_vars)
    s = length(mu) + n
    x_vars = all_indets[1:n]
    y_vars = all_indets[(n + 1):(n + m)]
    u_vars = all_indets[(n + m + 1):(n + m + u)]
    all_vars = all_indets[1:(n + m + u)]
    mu = all_indets[(n + m + u + 1):end] 
    
    Rjet = create_jet_ring(all_vars, mu, s+2)   
    gens_Rjet = gens(Rjet)
    z_aux = gens_Rjet[end-length(mu)] 
   
    x_eqs = [[add_to_var(unpack_fraction(xeq[1])[1], Rjet, 1), add_zero_to_vars(unpack_fraction(xeq[2])[1], Rjet) // add_zero_to_vars(unpack_fraction(xeq[2])[2], Rjet)] for xeq in x_eqs]
    y_eqs = [[add_to_var(unpack_fraction(yeq[1])[1], Rjet, 0), add_zero_to_vars(unpack_fraction(yeq[2])[1], Rjet) // add_zero_to_vars(unpack_fraction(yeq[2])[2], Rjet)] for yeq in y_eqs]

    eqs = vcat(x_eqs, y_eqs)
    Q = foldl(lcm,[unpack_fraction(ex[2])[2] for ex in eqs])
   
    not_int_cond_params = gens_Rjet[(end - length(mu) + 1):end]
    all_params = vcat(not_int_cond_params,gens_Rjet[1:n])   
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
    D1 = floor(BigInt, (length(params_to_assess) + 1) * 2 * d0 * s * (n + 1) * (1 + 2 * d0 * s) / (1 - p_local) )
    
    # (c, d) ---------------
    sample = sample_point(D1, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q)
    all_subs = sample[4]
    u_hat = sample[2]
    y_hat = sample[1]
    
    # (e) ------------------
    alpha = [1 for i in 1:n]
    beta = [0 for i in 1:m]
    Et = Array{fmpq_poly}(undef, 0)    
    x_theta_vars = all_params
    prolongation_possible = [1 for i in 1:m]
    
    # (f) ------------------
    all_x_theta_vars_subs = insert_zeros_to_vals(all_subs[1], all_subs[2])
    eqs_i_old = Array{fmpq_mpoly}(undef, 0)
    evl_old = Array{fmpq_mpoly}(undef, 0)
    while sum(prolongation_possible) > 0
        for i in 1:m
            if prolongation_possible[i] == 1 
                eqs_i = vcat(Et, Y[i][beta[i] + 1])
                evl     = [evaluate(eq, vcat(u_hat[1], y_hat[1]), vcat(u_hat[2],y_hat[2])) for eq in eqs_i if !(eq in eqs_i_old)]
                evl_old = vcat(evl_old, evl)
                JacX    = jacobi_matrix(evl_old, x_theta_vars, all_x_theta_vars_subs) 
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
    theta_l = Array{fmpq_mpoly}(undef,0)
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
        println("Not identifiable parameters:                      [", join(params_to_assess,", "), "]")
    else
        println("Locally identifiable parameters: [", join([get_order_var(th,non_jet_ring)[1] for th in theta_l],", "), "]")
        println("Not identifiable parameters:     [", join([get_order_var(th,non_jet_ring)[1] for th in setdiff(params_to_assess, theta_l)],", "), "]")

# 3. Randomize.
        println("Randomizing")
        # (a) ------------
        deg_variety =  foldl(*, [BigInt(total_degree(e)) for e in Et])
        D2 = floor(BigInt, 6 * length(theta_l) * deg_variety * (1 + 2 * d0 * maximum(beta)) / (1 - p_local) )
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
        vrs_sorted = vcat(sort([e for e in Et_x_vars], lt = (x, y) -> compare_diff_var(x, y, all_indets, n + m + u, s)), z_aux, sort(not_int_cond_params, rev=true))
# 4. Determine.
        println("GB computation")  
    
        Rjet_new, vrs_sorted = Singular.PolynomialRing(Singular.QQ, [string(v) for v in vrs_sorted], ordering=:degrevlex)
        theta_g = Array{spoly}(undef, 0)  
    
        Et_hat = [parent_ring_change(e, Rjet_new) for e in Et_hat]
        gb = GroebnerBasis.f4(Ideal(Rjet_new, vcat(Et_hat, parent_ring_change(z_aux * Q_hat, Rjet_new) - 1)), nthrds=64)
        println("Remainder computation")
        theta_l_new = [parent_ring_change(th, Rjet_new) for th in theta_l]
    
        theta = [parent_ring_change(th, Rjet_new) for th in params_to_assess]
        for i in 1:length(theta_l)
            if Singular.reduce(theta_l_new[i], gb) == Rjet_new(theta_hat[2][findfirst(isequal(theta_l[i]), theta_hat[1])]) 
                theta_g = vcat(theta_g, theta_l_new[i])
            end
        end
        println("\n=== Summary ===")
        println("Globally identifiable parameters:                 [", join([get_order_var(th,non_jet_ring)[1] for th in theta_g],", "), "]")
        println("Locally but not globally identifiable parameters: [", join([get_order_var(th,non_jet_ring)[1] for th in setdiff(theta_l_new,theta_g)],", "), "]")
        println("Not identifiable parameters:                      [", join([get_order_var(th,non_jet_ring)[1] for th in setdiff(theta, theta_l_new)],", "), "]")
        println("===============")
    end
end


