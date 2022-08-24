# ------------------------------------------------------------------------------
"""
    func sample_point(bound, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q)

Sample random values for parameters of the polynomial system.
"""
function sample_point(bound, x_vars, y_vars, u_variables, all_params, X_eq, Y_eq, Q, known_values=[], known_states_jet_form=[])
    local u_hat, theta_hat, all_subs

    s = length(all_params)
    y_hat_vars = Array{fmpq_mpoly}(undef, 0)
    y_hat_vals = Array{fmpq}(undef, 0)

    while true
        theta_hat = [fmpq(rnd) for rnd in rand(0:bound, s)]
        u_hat = [fmpq(rnd) for rnd in rand(0:bound, length(u_variables))]
        all_subs = [vcat(all_params, u_variables), vcat(theta_hat, u_hat)]
        for idx in 1:length(all_subs)
            if all_subs[idx][1] in known_states_jet_form
                all_subs[idx][2] = known_values[indexin(all_subs[idx][2], known_states_jet_form)]
            end
        end
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
