function get_weights(ode::SIAN.ODE{P}) where {P<:Nemo.MPolyElem{Nemo.fmpq}
    eqs, Q, x_eqs, y_eqs, x_vars, y_vars, u_vars, mu, gens_Rjet = get_equations(ode)

    n = length(x_vars)
    m = length(y_vars)
    u = length(u_vars)
    s = length(mu) + n

    non_jet_ring = ode.poly_ring
    all_indets = gens(non_jet_ring)
    z_aux = gens_Rjet[end-length(mu)]
    Rjet = gens_Rjet[1].parent

    n = length(x_vars)
    m = length(y_vars)
    u = length(u_vars)
    s = length(mu) + n

    current_level = 0
    visible_states = Dict{Int,Set{fmpq_mpoly}}(0=>Set{fmpq_mpoly}())
    for eq in values(ode.y_equations)
        visible_states[0] = union(visible_states[0], vars(eq))
    end

    differentiate_ = [true for each in ode.y_equations]

    for i in 1:s+1
        poly_d = SIAN.unpack_fraction(y_eqs[i][1] - y_eqs[i][2])[1]
        poly_d = SIAN.differentiate_all(poly_d, gens_Rjet, n + m + u, j - 1)
        leader = gens_Rjet[i+n+(n+m+u)*j]
        separant = SIAN.derivative(poly_d, leader)
        # poly_d := simplify(leader - subs(x_eqs, -(poly_d - separant * leader) / separant)):
    end

end