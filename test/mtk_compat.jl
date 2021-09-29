@testset "Check identifiability of `ODESystem` object" begin
    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t) [output = true]
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]

    de = ODESystem(eqs, t, name=:Test)
    correct = Dict{Any,Any}(
        "nonidentifiable" => [a12, a21, substitute(x1, t=>0), a01], 
        "locally_not_globally" => Num[], 
        "globally" => [substitute(x0, t=>0)]
    )
    output = identifiability_ode(de)
    @test isequal(correct["nonidentifiable"], output["nonidentifiable"])
    @test isequal(correct["locally_not_globally"], output["locally_not_globally"])
    @test isequal(correct["globally"], output["globally"])

    # --------------------------------------------------------------------------
        # # --------------------------------------------------------------------------
        # # ----------

    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t) [output = true]
    D = Differential(t)
    using SpecialFunctions

    eqs = [
        D(x0) ~ -(a01 + a21) * SpecialFunctions.erfc(x0) + a12 * x1,
        D(x1) ~ a21 * x0 - a12 * x1,
        y1 ~ x0
    ]

    de = ODESystem(eqs, t, name=:Test)
    
    @test_throws ArgumentError identifiability_ode(de )
    # ----------
    # @parameters a b c 
    # @variables t x1(t) x2(t) y(t) [output = true]
    # D = Differential(t)

    # eqs = [D(x1) ~ -a * x1 + x2 * b / (x1 + b / (c^2 - x2)), D(x2) ~ x2 * c^2 + x1, y ~ x2]
    # de = ODESystem(eqs, t, name=:Test)
    # correct = Dict(a => :globally, b => :globally, c => :locally)
    # to_check = [a, b, c]
    # @test isequal(correct, assess_identifiability(de, to_check))
end
