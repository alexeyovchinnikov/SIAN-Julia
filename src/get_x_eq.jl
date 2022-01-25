function get_x_eq(x_eqs::Vector{Vector{Nemo.AbstractAlgebra.RingElem}}, y_eqs::Vector{Vector{Nemo.AbstractAlgebra.RingElem}}, n::Int, m::Int, u::Int, gens_Rjet)
    X = Array{Nemo.fmpq_poly}(undef, 0)
    X_eq = Array{Nemo.fmpq_poly}(undef, 0)
    for i in 1:n
        X = vcat(X, [Array{Nemo.fmpq_poly}(undef, 0)])
        poly_d = SIAN.unpack_fraction(x_eqs[i][1] - x_eqs[i][2])[1]
        for j in 0:s+1
            if j > 0
                poly_d = SIAN.differentiate_all(poly_d, gens_Rjet, n + m + u, j)
            end
            leader = gens_Rjet[i+(n+m+u)*(j+1)]
            separant = SIAN.derivative(poly_d, leader)
            X[i] = vcat(X[i], poly_d)
            X_eq = vcat(X_eq, [[leader, -(poly_d - separant * leader) // separant]])
        end
    end
    return X, X_eq
end