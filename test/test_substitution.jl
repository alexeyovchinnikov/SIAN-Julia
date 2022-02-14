@testset "Testing weight substitution" begin
    R, (x, y, z, a, b) = SIAN.Nemo.PolynomialRing(SIAN.Nemo.QQ, [:x, :y, :z, :a, :b])

    Et_hat = [x + y, a * x + y * z, b - x^2 + y^3 * z]
    Et_hat_sub = [x^3 + y, a^2 * x^2 + y * z, b^2 - x^6 + y^3 * z]

    weights = Dict(a => 2, b => 2, x => 3)

    for i in 1:length(Et_hat)
        for _var in SIAN.Nemo.vars(Et_hat[i])
            Et_hat[i] = SIAN.StrucrturalIdentifiability.make_substitution(Et_hat[i], _var, _var^get(weights, _var, 1), parent(_var)(1))
        end
    end

    @test all(Et_hat_sub[i] == Et_hat[i] for i ∈ 1:length(Et_hat))

end