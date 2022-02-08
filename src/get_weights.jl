function get_weights(ode::SIAN.ODE{P}) where {P<:Nemo.MPolyElem{Nemo.fmpq}
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
    visible_states = Dict{Int,Set{fmpq_mpoly}}(0=>Set{fmpq_mpoly}())
    for eq in values(ode.y_equations)
        visible_states[0] = union(visible_states[0], vars(eq))
    end

    differentiate_ = [true for each in ode.y_equations]
    Y_eqs = Dict()
    vars_ = [a[1] for a in x_eqs]
    vals_ = [a[2] for a in x_eqs]
    for j in 1:s+1
        for i in 1:m
            poly_d = SIAN.unpack_fraction(y_eqs[i][1] - y_eqs[i][2])[1]
            poly_d = SIAN.differentiate_all(poly_d, gens_Rjet, n + m + u, j)
            leader = gens_Rjet[i+n+(n+m+u)*j]
            separant = SIAN.derivative(poly_d, leader)
            Y_eqs[leader] = -(poly_d - separant * leader) // separant

            candidates = select(x-> not (GetOrderVar(x)[1]  in y_vars), indets(y_eqs[i])):
        end
    end

end