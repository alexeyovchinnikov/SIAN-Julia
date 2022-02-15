@testset "Testing weight substitution" begin
    using StructuralIdentifiability: make_substitution
    @info "Case 1"
    R, (x, y, z, a, b) = SIAN.Nemo.PolynomialRing(SIAN.Nemo.QQ, [:x, :y, :z, :a, :b])
    Et_hat = [x + y, a * x + y * z, b - x^2 + y^3 * z]
    Et_hat_sub = [x^3 + y, a^2 * x^3 + y * z, b^2 - x^6 + y^3 * z]
    weights = Dict(a => 2, b => 2, x => 3)
    for i in 1:length(Et_hat)
        for _var in Set(SIAN.Nemo.vars(Et_hat[i]))
            println(_var, Et_hat[i])
            Et_hat[i] = make_substitution(Et_hat[i], _var, _var^get(weights, _var, 1), parent(_var)(1))
            println(_var, Et_hat[i])
        end
    end
    @test all(Et_hat_sub[i] == Et_hat[i] for i ∈ 1:length(Et_hat))

    @info "Case 2"
    R, (x, y, z, a, b, c) = SIAN.Nemo.PolynomialRing(SIAN.Nemo.QQ, [:x, :y, :z, :a, :b, :c])
    Et_hat = [x * y - a, 3 * x + y, 2 * x * y * z * a - z, x * y^2 + y^3 * z^4 - c^2 * b * x^2]
    Et_hat_sub = [x * y - a^2, 3 * x + y, 2 * x * y * z * a - z^4, x * y^2 + y^3 * z^20 - c^2 * b * x^2]
    weights = Dict(x => 1, a => 2, z => 5)
    for i in 1:length(Et_hat)
        for _var in Set(SIAN.Nemo.vars(Et_hat[i]))
            println(_var, Et_hat[i])
            Et_hat[i] = make_substitution(Et_hat[i], _var, _var^get(weights, _var, 1), parent(_var)(1))
            println(_var, Et_hat[i])
        end
    end

    @test all(Et_hat_sub[i] == Et_hat[i] for i ∈ 1:length(Et_hat))

    @info "Case 3"

    R, (x, y, z) = SIAN.Nemo.PolynomialRing(SIAN.Nemo.QQ, [:x, :y, :z, :a, :b])
    Et_hat = [x - y^3, x^3 + y, z - y^2, x + z^2]
    Et_hat_sub = [x - y^3, x + y, z - y^2, x + z^2]
    weights = Dict(x => 2, y => 2, z => 3)
    for i in 1:length(Et_hat)
        for _var in Set(SIAN.Nemo.vars(Et_hat[i]))
            println(_var, Et_hat[i])
            Et_hat[i] = make_substitution(Et_hat[i], _var, _var^get(weights, _var, 1), parent(_var)(1))
            println(_var, Et_hat[i])
        end
    end
    @test all(Et_hat_sub[i] == Et_hat[i] for i ∈ 1:length(Et_hat))

end