function get_weights(ode, non_identifiable_parameters)
    eqs, Q, x_eqs, y_eqs, x_vars, y_vars, u_vars, mu, all_indets, gens_Rjet = SIAN.get_equations(ode)

    non_jet_ring = ode.poly_ring
    all_indets = gens(non_jet_ring)
    z_aux = gens_Rjet[end-length(mu)]
    Rjet = gens_Rjet[1].parent

    n = length(x_vars)
    m = length(y_vars)
    u = length(u_vars)
    s = length(mu) + n

    current_level = 0
    visible_states = Dict{Int,Set{fmpq_mpoly}}(current_level => Set{fmpq_mpoly}())
    for eq in y_eqs
        numer, denom = SIAN.unpack_fraction(eq[2])
        # visible states must be in non-jet representation!
        visible_states[current_level] = union(visible_states[current_level], Set(SIAN.get_order_var(vn, non_jet_ring)[1] for vn in vars(numer)))
        visible_states[current_level] = union(visible_states[current_level], Set(SIAN.get_order_var(vd, non_jet_ring)[1] for vd in vars(denom)))
    end

    differentiate_ = [true for each in ode.y_equations]
    Y_eqs = Dict()
    X_eqs = Dict(a[1] => a[2] for a in x_eqs)
    seen_so_far = Set()
    for j in 1:s+1
        current_level += 1
        for i in 1:m
            if differentiate_[i] # if we need to differentiate current output
                # get the polynomial, Nemo takes care of common denominator
                poly_d = SIAN.unpack_fraction(y_eqs[i][1] - y_eqs[i][2])[1]

                # differentiate
                poly_d = SIAN.differentiate_all(poly_d, gens_Rjet, n + m + u, j)

                # get the leader (y[i]_k, the symbol for derivative of y)
                leader = gens_Rjet[i+n+(n+m+u)*j]

                # differentiate the polynomial w.r.t. to leader to get the "coefficient" of leader
                separant = SIAN.derivative(poly_d, leader)

                # write the equation
                Y_eqs[leader] = -(poly_d - separant * leader) // separant

                # split the numerator, the denominator
                numer, denom = SIAN.unpack_fraction(Y_eqs[leader])

                # substitute any derivatives of states with the original odes
                for v in vars(numer)
                    v_non_jet, v_order = SIAN.get_order_var(v, non_jet_ring)
                    if v_non_jet in x_vars && v_order > 0
                        numer = make_substitution(numer, v, SIAN.unpack_fraction(X_eqs[v])[1], SIAN.unpack_fraction(X_eqs[v])[2])
                        denom = make_substitution(denom, v, SIAN.unpack_fraction(X_eqs[v])[1], SIAN.unpack_fraction(X_eqs[v])[2])
                    end
                end

                # find states that newly appeared
                visible_states[current_level] = union(get(visible_states, current_level, Set{fmpq_mpoly}()), Set(SIAN.get_order_var(vn, non_jet_ring)[1] for vn in vars(numer)))
                visible_states[current_level] = union(get(visible_states, current_level, Set{fmpq_mpoly}()), Set(SIAN.get_order_var(vd, non_jet_ring)[1] for vd in vars(denom)))

                # add previous level to "what we have seen so far"-set
                union!(seen_so_far, visible_states[current_level-1])
                # remove states that were seen on prevous level
                visible_states[current_level] = setdiff(visible_states[current_level], seen_so_far)
                Y_eqs[leader] = numer // denom

                if length(visible_states[current_level]) == 0
                    differentiate_[i] = false
                end
                y_eqs[i] = [leader, Y_eqs[leader]]
            end
        end
        if sum(differentiate_) == 0
            break
        end
    end
    weights = Dict{fmpq_mpoly,Int64}()
    max_level = current_level - 1
    # TODO: filter non-identifiable parameters
    for (level, states) in visible_states
        for st in states
            if st in x_vars
                weights[st] = level + 1
            end
            if st in mu && level == max_level && !(st in non_identifiable_parameters)
                weights[st] = level + 1
            end
        end
    end
    return weights
end
