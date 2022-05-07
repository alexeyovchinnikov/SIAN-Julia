module SIAN

using Nemo
using StructuralIdentifiability: PreprocessODE, eval_at_nemo
using LinearAlgebra
using Singular
using GroebnerBasis
using MacroTools
using OrderedCollections
using ModelingToolkit
using Logging
include("util.jl")
include("ODE.jl")
include("max_poly_system.jl")
include("get_x_eq.jl")
include("get_y_eq.jl")
include("sample_point.jl")

export identifiability_ode, PreprocessODE
export @ODEmodel
export ODE
export Nemo, OrderedDict, Generic, macroexpand, macrohelper_extract_vars, macrohelper_clean, fmpq_mpoly, get_parameters

<<<<<<< HEAD
=======
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
    u_hat = [fmpq(rnd) for rnd in rand(0:bound, length(u_variables))]
    all_subs = [vcat(all_params, u_variables), vcat(theta_hat, u_hat)]
    if evaluate(Q, all_subs[1], all_subs[2]) == 0
      next
    end
    vls = insert_zeros_to_vals(all_subs[1], all_subs[2])
    for i in 0:(s+1)
      for j in 1:length(y_vars)
        eq = Y_eq[(j-1)*(s+2)+i+1][2]
        vl = evaluate(unpack_fraction(eq)[1], vls) // evaluate(unpack_fraction(eq)[2], vls)
        y_hat_vars = vcat(y_hat_vars, Y_eq[(j-1)*(s+2)+i+1][1])
        y_hat_vals = vcat(y_hat_vals, vl)
        vls[var_index(Y_eq[(j-1)*(s+2)+i+1][1])] = vl
        all_subs = [vcat(all_subs[1], Y_eq[(j-1)*(s+2)+i+1][1]), vcat(all_subs[2], vl)]
      end
      for j in 1:length(x_vars)
        eq = X_eq[(j-1)*(s+2)+i+1][2]
        vl = evaluate(unpack_fraction(eq)[1], vls) // evaluate(unpack_fraction(eq)[2], vls)
        all_subs = [vcat(all_subs[1], X_eq[(j-1)*(s+2)+i+1][1]), vcat(all_subs[2], vl)]
        vls[var_index(X_eq[(j-1)*(s+2)+i+1][1])] = vl
      end
    end
    break
  end
  return ([[y_hat_vars, y_hat_vals], [u_variables, u_hat], [all_params, theta_hat], all_subs])
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
function get_parameters(ode; initial_conditions=true)
  if initial_conditions
    return vcat(ode.parameters, ode.x_vars)
  else
    return ode.parameters
  end
end
>>>>>>> main

# ------------------------------------------------------------------------------
# Main Code
# ------------------------------------------------------------------------------
#

<<<<<<< HEAD
"""
    func identifiability_ode(ode, params_to_assess; p=0.99, p_mod=0, nthrds=1)

Perform identifiability check for a given `ode` system with respect to parameters in `params_to_assess` list.
    
## Input

  - `ode` - an ODE system returned by the `@ODEmodel` macro.
  - `params_to_assess` - an array of parameters returned by `get_parameters` function.
  - `p` - probability of correctness, default `0.99`.
  - `p_mod` - a prime characteristic, default `0`.
  - `nthrds` - number of threads for concurrency, default `1`.
  - `infolevel` - an integer, controls the verbosity of Groebner Basis computation, default `0` (no output).
                 See GroebnerBasis.jl documentation for details.
"""

function identifiability_ode(ode, params_to_assess; p = 0.99, p_mod = 0, nthrds = 1, infolevel = 0)

    println("Solving the problem")

    if p_mod != 0
        @warn "Using `p_mod!=0` does not guarantee the same probability of correctness but allows to run the program faster. This warning was raised by `p_mod = $p_mod`."
    end
    # 1.Construct the maximal system

    # (a) ---------------

    println("Constructing the maximal system")

    eqs, Q, x_eqs, y_eqs, x_vars, y_vars, u_vars, mu, all_indets, gens_Rjet = SIAN.get_equations(ode)

    non_jet_ring = ode.poly_ring
    z_aux = gens_Rjet[end-length(mu)]
    Rjet = gens_Rjet[1].parent

    n = length(x_vars)
    m = length(y_vars)
    u = length(u_vars)
    s = length(mu) + n

    not_int_cond_params = gens_Rjet[(end-length(ode.parameters)+1):end]
    all_params = vcat(not_int_cond_params, gens_Rjet[1:n])
    x_variables = gens_Rjet[1:n]
    for i in 1:(s+1)
        x_variables = vcat(x_variables, gens_Rjet[(i*(n+m+u)+1):(i*(n+m+u)+n)])
    end
    u_variables = gens_Rjet[(n+m+1):(n+m+u)]
    for i in 1:(s+1)
        u_variables = vcat(u_variables, gens_Rjet[((n+m+u)*i+n+m+1):((n+m+u)*(i+1))])
    end

    # (b,c) -------------
    X, X_eq = SIAN.get_x_eq(x_eqs, y_eqs, n, m, s, u, gens_Rjet)

    # (d,e) -----------
    Y, Y_eq = SIAN.get_y_eq(x_eqs, y_eqs, n, m, s, u, gens_Rjet)

    # 2.Truncate
    println("Truncating")

    # (a) -----------------------
    d0 = BigInt(maximum(vcat([total_degree(SIAN.unpack_fraction(Q * eq[2])[1]) for eq in eqs], total_degree(Q))))

    # (b) -----------------------  
    D1 = floor(BigInt, (length(params_to_assess) + 1) * 2 * d0 * s * (n + 1) * (1 + 2 * d0 * s) / (1 - p))

    # (c, d) ---------------
    sample = SIAN.sample_point(D1, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q)
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
    all_x_theta_vars_subs = SIAN.insert_zeros_to_vals(all_subs[1], all_subs[2])
    eqs_i_old = Array{fmpq_mpoly}(undef, 0)
    evl_old = Array{fmpq_mpoly}(undef, 0)
    while sum(prolongation_possible) > 0
        for i in 1:m
            if prolongation_possible[i] == 1
                eqs_i = vcat(Et, Y[i][beta[i]+1])
                evl = [evaluate(eq, vcat(u_hat[1], y_hat[1]), vcat(u_hat[2], y_hat[2])) for eq in eqs_i if !(eq in eqs_i_old)]
                evl_old = vcat(evl_old, evl)
                JacX = SIAN.jacobi_matrix(evl_old, x_theta_vars, all_x_theta_vars_subs)
                eqs_i_old = eqs_i
                if LinearAlgebra.rank(JacX) == length(eqs_i)
                    Et = vcat(Et, Y[i][beta[i]+1])
                    beta[i] = beta[i] + 1
                    # adding necessary X-equations
                    polys_to_process = vcat(Et, [Y[k][beta[k]+1] for k in 1:m])
                    while length(polys_to_process) != 0
                        new_to_process = Array{fmpq_mpoly}(undef, 0)
                        vrs = Set{fmpq_mpoly}()
                        for poly in polys_to_process
                            vrs = union(vrs, [v for v in vars(poly) if v in x_variables])
                        end
                        vars_to_add = Set{fmpq_mpoly}(v for v in vrs if !(v in x_theta_vars))
                        for v in vars_to_add
                            x_theta_vars = vcat(x_theta_vars, v)
                            ord_var = SIAN.get_order_var2(v, all_indets, n + m + u, s)
                            var_idx = var_index(ord_var[1])
                            poly = X[var_idx][ord_var[2]]
                            Et = vcat(Et, poly)
                            new_to_process = vcat(new_to_process, poly)
                            alpha[var_idx] = max(alpha[var_idx], ord_var[2] + 1)
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
        for j in (beta[i]+1):length(Y[i])
            to_add = true
            for v in SIAN.get_vars(Y[i][j], x_vars, all_indets, n + m + u, s)
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
    params_to_assess_ = [SIAN.add_to_var(param, Rjet, 0) for param in params_to_assess]
    Et_eval_base = [evaluate(e, vcat(u_hat[1], y_hat[1]), vcat(u_hat[2], y_hat[2])) for e in Et]
    for param_0 in params_to_assess_
        other_params = [v for v in x_theta_vars if v != param_0]
        Et_subs = [evaluate(e, [param_0], [evaluate(param_0, all_x_theta_vars_subs)]) for e in Et_eval_base]
        JacX = SIAN.jacobi_matrix(Et_subs, other_params, all_x_theta_vars_subs)
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
        println("Locally identifiable parameters: [", join([SIAN.get_order_var(th, non_jet_ring)[1] for th in theta_l], ", "), "]")
        println("Not identifiable parameters:     [", join([SIAN.get_order_var(th, non_jet_ring)[1] for th in setdiff(params_to_assess_, theta_l)], ", "), "]")

        # 3. Randomize.
        println("Randomizing")
        # (a) ------------
        deg_variety = foldl(*, [BigInt(total_degree(e)) for e in Et])
        D2 = floor(BigInt, 6 * length(theta_l) * deg_variety * (1 + 2 * d0 * maximum(beta)) / (1 - p))
        # (b, c) ---------
        sample = SIAN.sample_point(D2, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q)
        y_hat = sample[1]
        u_hat = sample[2]
        theta_hat = sample[3]

        # (d) ------------
        Et_hat = [evaluate(e, vcat(y_hat[1], u_hat[1]), vcat(y_hat[2], u_hat[2])) for e in Et]
        Et_x_vars = Set{fmpq_mpoly}()
        for poly in Et_hat
            Et_x_vars = union(Et_x_vars, Set(vars(poly)))
        end
        Q_hat = evaluate(Q, u_hat[1], u_hat[2])
        vrs_sorted = vcat(sort([e for e in Et_x_vars], lt = (x, y) -> SIAN.compare_diff_var(x, y, all_indets, n + m + u, s)), z_aux, sort(not_int_cond_params, rev = true))

        # 4. Determine.
        println("GB computation")

        if p_mod > 0
            Et_hat = [SIAN._reduce_poly_mod_p(e, p_mod) for e in Et_hat]
            z_aux = SIAN._reduce_poly_mod_p(z_aux, p_mod)
            Q_hat = SIAN._reduce_poly_mod_p(Q_hat, p_mod)
            Rjet_new, vrs_sorted = Singular.PolynomialRing(Singular.Fp(p_mod), [string(v) for v in vrs_sorted], ordering = :degrevlex)
        else
            Rjet_new, vrs_sorted = Singular.PolynomialRing(Singular.QQ, [string(v) for v in vrs_sorted], ordering = :degrevlex)
        end

        theta_g = Array{spoly}(undef, 0)
        Et_hat = [SIAN.parent_ring_change(e, Rjet_new) for e in Et_hat]
        gb = GroebnerBasis.f4(Ideal(Rjet_new, vcat(Et_hat, SIAN.parent_ring_change(z_aux * Q_hat, Rjet_new) - 1)), nthrds = nthrds, infolevel = infolevel)
        println("Remainder computation")

        if p_mod > 0
            theta_l_new = [SIAN.parent_ring_change(SIAN._reduce_poly_mod_p(th, p_mod), Rjet_new) for th in theta_l]

            for i in 1:length(theta_l)
                if Singular.reduce(theta_l_new[i], gb) == SIAN.parent_ring_change(SIAN._reduce_poly_mod_p(Rjet(theta_hat[2][findfirst(isequal(theta_l[i]), theta_hat[1])]), p_mod), Rjet_new)
                    theta_g = vcat(theta_g, theta_l_new[i])
                end
            end
        else
            theta_l_new = [SIAN.parent_ring_change(th, Rjet_new) for th in theta_l]

            for i in 1:length(theta_l)
                if Singular.reduce(theta_l_new[i], gb) == SIAN.parent_ring_change(Rjet(theta_hat[2][findfirst(isequal(theta_l[i]), theta_hat[1])]), Rjet_new)
                    theta_g = vcat(theta_g, theta_l_new[i])
                end
            end

        end
        result = Dict(
            "globally" => Set(SIAN.get_order_var(th, non_jet_ring)[1] for th in theta_g),
            "locally_not_globally" => Set(SIAN.get_order_var(th, non_jet_ring)[1] for th in setdiff(theta_l_new, theta_g)),
            "nonidentifiable" => Set(SIAN.get_order_var(th, non_jet_ring)[1] for th in setdiff(params_to_assess_, theta_l))
        )
        println("\n=== Summary ===")
        println("Globally identifiable parameters:                 [", join(result["globally"], ", "), "]")
        println("Locally but not globally identifiable parameters: [", join(result["locally_not_globally"], ", "), "]")
        println("Not identifiable parameters:                      [", join(result["nonidentifiable"], ", "), "]")
        println("===============")
        return result
    end
end


function identifiability_ode(ode::ModelingToolkit.ODESystem, params_to_assess = []; p = 0.99, p_mod = 0, nthrds = 1)
    if any(ModelingToolkit.isoutput(eq.lhs) for eq in ModelingToolkit.equations(ode))
        # @info "Measured quantities are not provided, trying to find the outputs in input ODE."
        measured_quantities = filter(eq -> (ModelingToolkit.isoutput(eq.lhs)), ModelingToolkit.equations(ode))
    else
        throw(error("Measured quantities (output functions) were not provided and no outputs were found."))
    end
    ode_prep, input_syms, gens_ = PreprocessODE(ode, measured_quantities)
    t = ModelingToolkit.arguments(ModelingToolkit.states(ode)[1])[1]
    if length(params_to_assess) == 0
        params_to_assess_ = SIAN.get_parameters(ode_prep)
        nemo2mtk = Dict(gens_ .=> input_syms)
    else
        params_to_assess_ = [eval_at_nemo(each, Dict(syms .=> gens_)) for each in params_to_assess]
        nemo2mtk = Dict(params_to_assess_ .=> params_to_assess)
    end
=======
function identifiability_ode(ode::ModelingToolkit.ODESystem, params_to_assess=[]; p=0.99, p_mod=0, nthrds=1)
  if any(ModelingToolkit.isoutput(eq.lhs) for eq in ModelingToolkit.equations(ode))
    # @info "Measured quantities are not provided, trying to find the outputs in input ODE."
    measured_quantities = filter(eq -> (ModelingToolkit.isoutput(eq.lhs)), ModelingToolkit.equations(ode))
  else
    throw(error("Measured quantities (output functions) were not provided and no outputs were found."))
  end
  ode_prep, input_syms, gens_ = PreprocessODE(ode, measured_quantities)
  t = ModelingToolkit.arguments(ModelingToolkit.states(ode)[1])[1]
  if length(params_to_assess) == 0
    params_to_assess_ = SIAN.get_parameters(ode_prep)
    nemo2mtk = Dict(gens_ .=> input_syms)
  else
    params_to_assess_ = [eval_at_nemo(each, Dict(syms .=> gens_)) for each in params_to_assess]
    nemo2mtk = Dict(params_to_assess_ .=> params_to_assess)
  end
>>>>>>> main

  res = identifiability_ode(ode_prep, params_to_assess_; p=p, p_mod=p_mod, nthrds=1)

  @info "Post-Processing: Converting Nemo output to ModelingToolkit types"
  out = Dict()
  for (id_type, pars) in pairs(res)
    out[id_type] = [ModelingToolkit.Num(substitute(nemo2mtk[each], t => 0)) for each in pars]
  end
  return out
end


<<<<<<< HEAD
=======
Perform identifiability check for a given `ode` system with respect to parameters in `params_to_assess` list.
    
## Input

  - `ode` - an ODE system returned by the `@ODEmodel` macro.
  - `params_to_assess` - an array of parameters returned by `get_parameters` function.
  - `p` - probability of correctness, default `0.99`.
  - `p_mod` - a prime characteristic, default `0`.
  - `nthrds` - number of threads for concurrency, default `1`.
  - `infolevel` - an integer, controls the verbosity of Groebner Basis computation, default `0` (no output).
                 See GroebnerBasis.jl documentation for details.
"""
function identifiability_ode(ode, params_to_assess; p=0.99, p_mod=0, nthrds=1, infolevel=0)

  @info "Solving the problem"

  if p_mod != 0
    @warn "Using `p_mod!=0` does not guarantee the same probability of correctness but allows to run the program faster. This warning was raised by `p_mod = $p_mod`."
  end
  # 1.Construct the maximal system

  # (a) ---------------

  @info "Constructing the maximal system"

  non_jet_ring = ode.poly_ring
  all_indets = gens(non_jet_ring)
  x_vars = ode.x_vars
  y_vars = ode.y_vars
  u_vars = ode.u_vars
  mu = ode.parameters

  n = length(x_vars)
  m = length(y_vars)
  u = length(u_vars)
  s = length(mu) + n

  Rjet = create_jet_ring(vcat(x_vars, y_vars, u_vars), mu, s + 2)
  gens_Rjet = gens(Rjet)
  z_aux = gens_Rjet[end-length(mu)]
  x_eqs = collect(values(ode.x_equations))
  y_eqs = collect(values(ode.y_equations))

  x_eqs = [[add_to_var(x_vars[i], Rjet, 1), add_zero_to_vars(unpack_fraction(x_eqs[i])[1], Rjet) // add_zero_to_vars(unpack_fraction(x_eqs[i])[2], Rjet)] for i in 1:n]
  y_eqs = [[add_to_var(y_vars[i], Rjet, 0), add_zero_to_vars(unpack_fraction(y_eqs[i])[1], Rjet) // add_zero_to_vars(unpack_fraction(y_eqs[i])[2], Rjet)] for i in 1:m]

  eqs = vcat(x_eqs, y_eqs)
  Q = foldl(lcm, [unpack_fraction(ex[2])[2] for ex in eqs])

  not_int_cond_params = gens_Rjet[(end-length(ode.parameters)+1):end]
  all_params = vcat(not_int_cond_params, gens_Rjet[1:n])
  x_variables = gens_Rjet[1:n]
  for i in 1:(s+1)
    x_variables = vcat(x_variables, gens_Rjet[(i*(n+m+u)+1):(i*(n+m+u)+n)])
  end
  u_variables = gens_Rjet[(n+m+1):(n+m+u)]
  for i in 1:(s+1)
    u_variables = vcat(u_variables, gens_Rjet[((n+m+u)*i+n+m+1):((n+m+u)*(i+1))])
  end

  # (b,c) -------------
  X = Array{fmpq_poly}(undef, 0)
  X_eq = Array{fmpq_poly}(undef, 0)
  for i in 1:n
    X = vcat(X, [Array{fmpq_poly}(undef, 0)])
    poly_d = unpack_fraction(x_eqs[i][1] - x_eqs[i][2])[1]
    for j in 0:s+1
      if j > 0
        poly_d = differentiate_all(poly_d, gens_Rjet, n + m + u, j)
      end
      leader = gens_Rjet[i+(n+m+u)*(j+1)]
      separant = derivative(poly_d, leader)
      X[i] = vcat(X[i], poly_d)
      X_eq = vcat(X_eq, [[leader, -(poly_d - separant * leader) // separant]])
    end
  end

  # (d,e) -----------
  Y = Array{fmpq_poly}(undef, 0)
  Y_eq = Array{fmpq_poly}(undef, 0)
  for i in 1:m
    Y = vcat(Y, [Array{fmpq_poly}(undef, 0)])
    poly_d = unpack_fraction(y_eqs[i][1] - y_eqs[i][2])[1]
    for j in 0:s+1
      if j > 0
        poly_d = differentiate_all(poly_d, gens_Rjet, n + m + u, j - 1)
      end
      leader = gens_Rjet[i+n+(n+m+u)*j]
      separant = derivative(poly_d, leader)
      Y[i] = vcat(Y[i], poly_d)
      Y_eq = vcat(Y_eq, [[leader, -(poly_d - separant * leader) // separant]])
    end
  end

  # 2.Truncate
  @info "Truncating"

  # (a) -----------------------
  d0 = BigInt(maximum(vcat([total_degree(unpack_fraction(Q * eq[2])[1]) for eq in eqs], total_degree(Q))))

  # (b) -----------------------  
  D1 = floor(BigInt, (length(params_to_assess) + 1) * 2 * d0 * s * (n + 1) * (1 + 2 * d0 * s) / (1 - p))

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
        eqs_i = vcat(Et, Y[i][beta[i]+1])
        evl = [evaluate(eq, vcat(u_hat[1], y_hat[1]), vcat(u_hat[2], y_hat[2])) for eq in eqs_i if !(eq in eqs_i_old)]
        evl_old = vcat(evl_old, evl)
        JacX = jacobi_matrix(evl_old, x_theta_vars, all_x_theta_vars_subs)
        eqs_i_old = eqs_i
        if LinearAlgebra.rank(JacX) == length(eqs_i)
          Et = vcat(Et, Y[i][beta[i]+1])
          beta[i] = beta[i] + 1
          # adding necessary X-equations
          polys_to_process = vcat(Et, [Y[k][beta[k]+1] for k in 1:m])
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
              poly = X[var_idx][ord_var[2]]
              Et = vcat(Et, poly)
              new_to_process = vcat(new_to_process, poly)
              alpha[var_idx] = max(alpha[var_idx], ord_var[2] + 1)
            end
            polys_to_process = new_to_process
          end
        else
          prolongation_possible[i] = 0
        end
      end
    end
  end

  @info "Assessing local identifiability"

  max_rank = length(Et)
  for i in 1:m
    for j in (beta[i]+1):length(Y[i])
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
    @info "=== Summary ==="
    @info "Globally identifiable parameters:                 []"
    @info "Locally but not globally identifiable parameters: []"
    @info "Not identifiable parameters:                      [$(join(params_to_assess, ", "))]"
    return Dict("locally_identifiable" => [], "globally_identifiable" => [], "non_identifiable" => Set(SIAN.get_order_var(th, non_jet_ring)[1] for th in params_to_assess))
  else
    @info "Locally identifiable parameters: [$(join([get_order_var(th, non_jet_ring)[1] for th in theta_l], ", "))]"
    @info "Not identifiable parameters:     [$(join([get_order_var(th, non_jet_ring)[1] for th in setdiff(params_to_assess, theta_l)], ", "))]"

    # 3. Randomize.
    @info "Randomizing"
    # (a) ------------
    deg_variety = foldl(*, [BigInt(total_degree(e)) for e in Et])
    D2 = floor(BigInt, 6 * length(theta_l) * deg_variety * (1 + 2 * d0 * maximum(beta)) / (1 - p))
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
    @info "GB computation"

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
    gb = GroebnerBasis.f4(Ideal(Rjet_new, vcat(Et_hat, parent_ring_change(z_aux * Q_hat, Rjet_new) - 1)), nthrds=nthrds, infolevel=infolevel)
    @info "Remainder computation"

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
    result = Dict(
      "globally" => Set([get_order_var(th, non_jet_ring)[1] for th in theta_g]),
      "locally_not_globally" => Set([get_order_var(th, non_jet_ring)[1] for th in setdiff(theta_l_new, theta_g)]),
      "nonidentifiable" => Set([get_order_var(th, non_jet_ring)[1] for th in setdiff(params_to_assess, theta_l)])
    )
    @info "=== Summary ==="
    @info "Globally identifiable parameters:                 [$(join(result["globally"], ", "))]"
    @info "Locally but not globally identifiable parameters: [$(join(result["locally_not_globally"], ", "))]"
    @info "Not identifiable parameters:                      [$(join(result["nonidentifiable"], ", "))]"
    @info "==============="
    return result
  end
end
>>>>>>> main
end
