function get_equations(ode::SIAN.ODE{P}) where {P<:Nemo.MPolyElem{Nemo.fmpq}}
    non_jet_ring = ode.poly_ring
    all_indets = Nemo.gens(non_jet_ring)
    x_vars = ode.x_vars
    y_vars = ode.y_vars
    u_vars = ode.u_vars
    mu = ode.parameters

    n = length(x_vars)
    m = length(y_vars)
    u = length(u_vars)
    s = length(mu) + n

    Rjet = SIAN.create_jet_ring(vcat(x_vars, y_vars, u_vars), mu, s + 2)
    gens_Rjet = Nemo.gens(Rjet)
    z_aux = gens_Rjet[end-length(mu)]

    x_eqs = Array{Array{Nemo.AbstractAlgebra.RingElem,1},1}(undef, n)
    y_eqs = Array{Array{Nemo.AbstractAlgebra.RingElem,1},1}(undef, m)
    for i in 1:n
        numer, denom = SIAN.unpack_fraction(ode.x_equations[x_vars[i]])
        x_eqs[i] = [SIAN.add_to_var(x_vars[i], Rjet, 1), SIAN.add_zero_to_vars(numer, Rjet) // SIAN.add_zero_to_vars(denom, Rjet)]
    end
    for i in 1:m
        numer, denom = SIAN.unpack_fraction(ode.y_equations[y_vars[i]])
        y_eqs[i] = [SIAN.add_to_var(y_vars[i], Rjet, 0), SIAN.add_zero_to_vars(numer, Rjet) // SIAN.add_zero_to_vars(denom, Rjet)]
    end

    eqs = vcat(x_eqs, y_eqs)
    Q = foldl(lcm, [SIAN.unpack_fraction(ex[2])[2] for ex in eqs])

    return eqs, Q, x_eqs, y_eqs, x_vars, y_vars, u_vars, mu, gens_Rjet
end
