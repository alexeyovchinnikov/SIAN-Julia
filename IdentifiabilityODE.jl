println("Loading packages")
using Oscar
using LinearAlgebra
using Singular
using GroebnerBasis

function str_to_var(s, ring)
     return gens(ring)[findfirst(v -> (string(v) == s), gens(ring))]
end

function MakeDerivative(var_name, der_order)
  return(string(var_name, "_",der_order))
end

function JetVars2D(var_list, max_ord)
  return [MakeDerivative(var,der_order) for var in var_list , der_order in 0:max_ord]
end

function AddOneToVar(vr, ring)
  return str_to_var(MakeDerivative(vr,1),ring)
end

function JetPolynomialRing(var_list, max_ord)
  local jet_vars
  jet_vars = vec(JetVars2D(var_list, max_ord)) 
  return PolynomialRing(QQ, jet_vars) 
end

function CreateJetRing(var_list, param_list, max_ord)
   varnames = vcat(vec(["$(s)_$i" for s in var_list, i in 0:max_ord]), "z_aux", ["$(s)_0" for s in param_list])
   return Oscar.PolynomialRing(Oscar.QQ, varnames)[1]
end

function AddZeroToVar(vr, ring)
  return str_to_var(MakeDerivative(vr,0),ring)
end

function DifferentiateOnce(diff_poly, var_list, max_ord)
  local result, i, j, jet_vars_2D 
  jet_vars_2D = JetVars2D(var_list, max_ord+1)
  result = 0
  for i in 1:length(var_list)
  for j in 1:max_ord+1
      result = result + derivative(diff_poly,str_to_var(jet_vars_2D[i, j], parent(diff_poly))) * str_to_var(jet_vars_2D[i, j+1], parent(diff_poly))
   end
  end
  return(result)
end

function DifferentiateMany(diff_poly, var_list, max_ords, ord)
  local result, i
  result = diff_poly
  for i in 1:ord 
    result = DifferentiateOnce(result, var_list, max_ords+i-1)
  end
  return(result)
end

function unpack_fraction(f)
  if applicable(numerator, f)
    return (numerator(f), denominator(f))
  end
  return (f, parent(f)(1))
end


function eval_frac_num(frac, vars, vals)
  if iszero(frac)
    return 0
  else
    fr = unpack_fraction(frac)
    vars = [unpack_fraction(vr)[1] for vr in vars]
    vals = [unpack_fraction(vl)[1] for vl in vals]
    return(fmpq(Oscar.lead(evaluate(fr[1],vars,vals))//Oscar.lead(evaluate(fr[2],vars,vals))))
  end
end

function eval_frac(frac, vars, vals)
  fr = unpack_fraction(frac)
  vars = [unpack_fraction(vr)[1] for vr in vars]
  vals = [unpack_fraction(vl)[1] for vl in vals]
  return(evaluate(fr[1],vars,vals)//evaluate(fr[2],vars,vals))
end

function is_numeric(frac)
  fr = unpack_fraction(frac)
  return(total_degree(fr[1])==0 & total_degree(fr[2])==0)
end

function SamplePoint(bound, x_vars, y_vars, u_vars, mu, X_eq, Y_eq, Q)
  local n, m, s, all_params, all_vars, theta_hat, u_variables, u_hat, x_hat, y_hat, v, poly, all_subs, to_compute
  n = length(x_vars)
  m = length(y_vars)
  s = length(mu) + n
  all_params = [AddZeroToVar(var,parent(Q)) for var in vcat(mu,x_vars)]
  all_vars   = vcat(x_vars, y_vars, u_vars)
  
  while true
     theta_hat = [parent(Q)(rnd) for rnd in rand(0:bound, length(all_params))] 
     u_variables = vec([str_to_var(MakeDerivative(u_vars[i], j),parent(Q)) for i in 1:length(u_vars),j in 0:(s+1)])
     u_hat =  [parent(Q)(rnd) for rnd in rand(0:bound, length(u_variables))]    

     all_subs = [vcat(all_params,u_variables), vcat(theta_hat, u_hat)]
     if evaluate(Q, all_subs[1], all_subs[2]) == 0
        next
     end
     to_compute = vcat(X_eq, Y_eq)
 
     while length(to_compute) != 0 
       to_compute = [[e[1],eval_frac(e[2], all_subs[1], all_subs[2])] for e in to_compute]
       new_subs = [[e[1],e[2]] for e in to_compute if is_numeric(e[2])]
       all_subs = [vcat(all_subs[1], [new_subs[i][1] for i in 1:length(new_subs)]), vcat(all_subs[2],[new_subs[i][2] for i in 1:length(new_subs)])]
       to_compute = [[e[1],e[2]] for e in to_compute if !is_numeric(e[2])]
     end

     y_hat = [[e[1] for e in Y_eq],[fmpq(eval_frac_num(e[2],all_subs[1],all_subs[2])) for e in Y_eq]]
     break
  end
  return([y_hat, [u_variables,u_hat], [all_params,theta_hat], all_subs])
end


function JacobiMatrix(pol_arr, vars, vals)
  m = Oscar.MatrixSpace(Oscar.QQ,length(pol_arr),length(vars))()
  for i in 1:length(pol_arr)
    for j in 1:length(vars)
      m[i,j] = eval_frac_num(derivative(pol_arr[i],vars[j]),vars,vals)
    end
  end
  return(m)
end

function GetOrderVar(diff_var,non_jet_ring)
  rex = match(r"^(.*_)([0-9]+)$", string(diff_var))
  if rex === nothing
    return(["", ""])
  else
     return([str_to_var(first(rex[1],length(rex[1])-1),non_jet_ring), parse(Int,rex[2])])
  end
end

function GetVars(diff_poly, var_list)
  local vrs
  vrs = union(Set{fmpq_mpoly}(vars(unpack_fraction(diff_poly)[1])), Set{fmpq_mpoly}(vars(unpack_fraction(diff_poly)[2])))
  return [v for v in vrs if GetOrderVar(v,parent(var_list[1]))[1] in var_list]
end

function CompareDiffVar(dvl, dvr, non_jet_ring)
  local vl, vr, hl, hr;
  vl, hl = GetOrderVar(dvl, non_jet_ring)
  vr, hr = GetOrderVar(dvr, non_jet_ring)
  if hl != hr 
    return (hl > hr)
  end
  if length(string(vl)) != length(string(vr)) 
    return (length(string(vl)) > length(string(vr)))
  end
  return (vr>=vl)
end 

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

function AddZeroToVars(poly::MPolyElem, new_ring::MPolyRing)
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

############# Main Code
function IdentifiabilityODE(x_eqs, y_eqs, u_vars, params_to_assess, int_cond_to_assess, p)
#equations will be arrays with [1] for lhs and [2] for rhs

   println("Solving the problem")
   # (a) ---------------
  non_jet_ring = parent(x_eqs[1][1])
  x_vars = [xv[1] for xv in x_eqs]
  y_vars = [yv[1] for yv in y_eqs]

# 1.Construct the maximal system
   println("Constructing the maximal system")
   # (a) ---------------
   params_to_assess = vcat(params_to_assess, int_cond_to_assess) 
   all_vars   = vcat(x_vars, y_vars, u_vars)
   #computing all symbols in the equations
   all_indets = vars(x_eqs[1][2])
   for i in 2:length(sigma_x)
      all_indets = union(all_indets,vars(sigma_x[i][2]))
   end
   for i in 1:length(sigma_y)
      all_indets = union(all_indets,vars(sigma_y[i][2]))
   end
   #finding which ones are the parameters
   mu = setdiff(all_indets, all_vars)
   p_local = p + length(params_to_assess) * 10^(-18)

   n = length(x_vars)
   m = length(y_vars)
   s = length(mu) + n

   Rjet = CreateJetRing(all_vars, mu, s+2)   

   z_aux = gens(Rjet)[end-length(mu)] 

   x_eqs = [[AddOneToVar(xeq[1],Rjet),AddZeroToVars(xeq[2],Rjet)] for xeq in x_eqs]
   y_eqs = [[AddZeroToVar(yeq[1],Rjet),AddZeroToVars(yeq[2],Rjet)] for yeq in y_eqs]


   eqs = vcat(x_eqs, y_eqs)
   Q = foldl(lcm,[unpack_fraction(ex[2])[2] for ex in eqs])

   all_params = [AddZeroToVar(var,Rjet) for var in vcat(mu,x_vars)]

   # (b,c) -------------
   X = Array{fmpq_poly}(undef,0)
   X_eq = Array{fmpq_poly}(undef,0)
   for i in 1:n
      X = vcat(X, [Array{fmpq_poly}(undef,0)])
      poly = unpack_fraction(x_eqs[i][1] - x_eqs[i][2])[1]  
      for j in 0:s+1
         poly_d = DifferentiateMany(poly, all_vars, 1, j) 
         leader = str_to_var(MakeDerivative(x_vars[i], j+1),parent(all_params[1]))
         separant = derivative(poly_d, leader)
         X[i] = vcat(X[i], poly_d)
         X_eq = vcat(X_eq, [[leader,-(poly_d - separant*leader)//separant]])
      end
   end
   
   # (d,e) -----------
   Y = Array{fmpq_poly}(undef,0)
   Y_eq = Array{fmpq_poly}(undef,0)
   for i in 1:m
      Y = vcat(Y, [Array{fmpq_poly}(undef,0)])
      poly = unpack_fraction(y_eqs[i][1] - y_eqs[i][2])[1]                
      for j in 0:s+1
         poly_d = DifferentiateMany(poly, all_vars, 0, j) 
         leader = str_to_var(MakeDerivative(y_vars[i], j),parent(all_params[1]))
         separant = derivative(poly_d, leader)
         Y[i] = vcat(Y[i], poly_d)
         Y_eq = vcat(Y_eq, [[leader,-(poly_d - separant*leader)//separant]])
      end
   end   
  
# 2.Truncate
   println("Truncating")
   # (a) -----------------------
   d0 = maximum(vcat([total_degree(unpack_fraction(Q*eq[2])[1]) for eq in eqs], total_degree(Q)))
   
   # (b) -----------------------  
   D1 = floor(BigInt, (length(params_to_assess) + 1) * 2 * d0 * s * (n + 1) * (1 + 2 * d0 * s) / (1 - p_local) )

   # (c, d) ---------------
   sample = SamplePoint(D1, x_vars, y_vars, u_vars, mu, X_eq, Y_eq, Q)
   all_subs = sample[4]
  u_hat = sample[2]
  y_hat = sample[1]
 
  # (e) ------------------
  alpha = [1 for i in 1:n]
  beta = [0 for i in 1:m]
  Et = Array{fmpq_poly}(undef,0)    
  x_theta_vars = all_params
  prolongation_possible = [1 for i in 1:m]

  # (f) ------------------
  while sum(prolongation_possible) > 0
    for i in 1:m 
      if prolongation_possible[i] == 1 
        eqs_i = vcat(Et, Y[i][beta[i] + 1])
        x_theta_vars_sub = [eval_frac(var,all_subs[1],all_subs[2]) for var in x_theta_vars]
        JacX = JacobiMatrix([eval_frac(eq, vcat(u_hat[1], y_hat[1]), vcat(u_hat[2],y_hat[2])) for eq in eqs_i], x_theta_vars, x_theta_vars_sub)
        if LinearAlgebra.rank(JacX) == length(eqs_i)
          Et = vcat(Et, Y[i][beta[i] + 1])
          beta[i] = beta[i] + 1
          # adding necessary X-equations
          polys_to_process = vcat(Et, [Y[k][beta[k] + 1] for k in 1:m])
          while length(polys_to_process) != 0 
            new_to_process = Array{fmpq_mpoly}(undef,0)
            vrs = Set{fmpq_mpoly}()
            for poly in polys_to_process
              vrs = union(vrs, GetVars(poly, x_vars))
            end
            vars_to_add = Set{fmpq_mpoly}(v for v in vrs if !(v in x_theta_vars)) 
            for v in vars_to_add
              x_theta_vars = vcat(x_theta_vars, v)
              ord_var = GetOrderVar(v,non_jet_ring)
              var_index = findfirst(x->x==ord_var[1], x_vars)
              poly = X[ var_index ][ ord_var[2] ]
              Et = vcat(Et, poly)
              new_to_process = vcat(new_to_process, poly)
              alpha[ var_index ] = max(alpha[ var_index ], ord_var[2] + 1)
            end
            polys_to_process = new_to_process
          end
        else
          prolongation_possible[i] = 0
        end
      end 
    end 
  end
  # is used for assessing local identifiabilty
  max_rank = length(Et)
  for i in 1:m
    for j in (beta[i] + 1):length(Y[i])
      to_add = true
      for v in GetVars(Y[i][j], x_vars)
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
  for param in params_to_assess
    param_0 = AddZeroToVar(param,Rjet)
    other_params = [v for v in x_theta_vars if v != param_0]
    Et_subs = [eval_frac(e, vcat(u_hat[1],y_hat[1],param_0), vcat(u_hat[2],y_hat[2],eval_frac(param_0,all_subs[1],all_subs[2]))) for e in Et]
    JacX = JacobiMatrix(Et_subs, other_params, [eval_frac(var,all_subs[1],all_subs[2]) for var in other_params])
    if LinearAlgebra.rank(JacX) != max_rank 
      theta_l = vcat(theta_l, param)
    end
  end
  if length(theta_l) > 0
    # 3. Randomize.
    println("Randomizing")
    # (a) ------------
    deg_variety = foldl(*, [total_degree(e) for e in Et])
    D2 = floor(BigInt, 6 * length(theta_l) * deg_variety * (1 + 2 * d0 * maximum(beta)) / (1 - p_local) )

    # (b, c) ---------
    sample = SamplePoint(D2, x_vars, y_vars, u_vars, mu, X_eq, Y_eq, Q)
    y_hat = sample[1]
    u_hat = sample[2]
    theta_hat = sample[3]

    # (d) ------------
    Et_hat = [evaluate(e, vcat(y_hat[1], u_hat[1]), vcat(y_hat[2], u_hat[2])) for e in Et]
    Et_x_vars = Set{fmpq_mpoly}()
    for poly in Et_hat
      Et_x_vars = union(Et_x_vars, GetVars(poly, x_vars))
    end
    Q_hat = evaluate(Q, u_hat[1], u_hat[2])

    vrs_sorted = vcat(sort([e for e in Et_x_vars], lt=(x,y)->CompareDiffVar(x,y,non_jet_ring)), z_aux, sort([AddZeroToVar(pr,Rjet) for pr in mu],rev=true))
  
    # 4. Determine.
  
    theta_l = [AddZeroToVar(th,Rjet) for th in theta_l]
    
    
    Rjet_new, vrs_sorted = Singular.PolynomialRing(Singular.QQ, [string(v) for v in vrs_sorted], ordering=:degrevlex)
    theta_g = Array{spoly}(undef,0)  
    Et_hat = [parent_ring_change(e, Rjet_new) for e in Et_hat]
    println("GB computation")
    gb = GroebnerBasis.f4(Ideal(Rjet_new,vcat(Et_hat, parent_ring_change(z_aux * Q_hat, Rjet_new) - 1)))

    println("Remainder computation")
    theta_l_new = [parent_ring_change(th, Rjet_new) for th in theta_l]
    
    theta = [AddZeroToVar(th, Rjet_new) for th in params_to_assess]
    for i in 1:length(theta_l)
      if Singular.reduce(theta_l_new[i], gb) == parent_ring_change(theta_hat[2][findfirst(isequal(theta_l[i]),theta_hat[1])], Rjet_new) 
        theta_g = vcat(theta_g, theta_l_new[i])
      end
    end

    println("\n=== Summary ===")
    println("Globally identifiable parameters:                 [", join([GetOrderVar(th,non_jet_ring)[1] for th in theta_g],","),"]")
    println("Locally but not globally identifiable parameters: [", join([GetOrderVar(th,non_jet_ring)[1] for th in setdiff(theta_l_new,theta_g)],","),"]")
    println("Not identifiable parameters:                      [", join([GetOrderVar(th,non_jet_ring)[1] for th in setdiff(theta, theta_l_new)],","),"]")
    println("===============")
  else
    println("\n=== Summary ===")
    println("Globally identifiable parameters:                 []")
    println("Locally but not globally identifiable parameters: []")
    println("Not identifiable parameters:                      [", join(params_to_assess,","),"]")
  end
end


