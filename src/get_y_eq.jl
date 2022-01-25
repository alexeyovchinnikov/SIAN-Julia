function get_y_eq(x_eqs::Vector{Vector{Nemo.AbstractAlgebra.RingElem}}, y_eqs::Vector{Vector{Nemo.AbstractAlgebra.RingElem}}, n::Int, m::Int, u::Int, gens_Rject)
    Y = Array{Nemo.fmpq_poly}(undef, 0)
    Y_eq = Array{Nemo.fmpq_poly}(undef, 0)
    for i in 1:m
        Y = vcat(Y, [Array{Nemo.fmpq_poly}(undef, 0)])
        poly_d = SIAN.unpack_fraction(y_eqs[i][1] - y_eqs[i][2])[1]
        for j in 0:s+1
            if j > 0
                poly_d = SIAN.differentiate_all(poly_d, gens_Rjet, n + m + u, j - 1)
            end
            leader = gens_Rjet[i+n+(n+m+u)*j]
            separant = SIAN.derivative(poly_d, leader)
            Y[i] = vcat(Y[i], poly_d)
            Y_eq = vcat(Y_eq, [[leader, -(poly_d - separant * leader) // separant]])
        end
    end
    return Y, Y_eq
end