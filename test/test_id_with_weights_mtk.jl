@testset "Testing identifiability with weights for MTK systems" begin
    @independent_variables t
    @parameters a01 a21 a12
    @variables x0(t) x1(t) y1(t) [output = true]
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]

    de = ODESystem(eqs, t, name = :Test)
    correct = Dict{Any,Any}(
        "nonidentifiable" => [a12, a21, substitute(x1, t => 0), a01],
        "locally_not_globally" => Num[],
        "globally" => [substitute(x0, t => 0)]
    )
    output = identifiability_ode(de; weighted_ordering = true)

    output_no_weight = identifiability_ode(de; weighted_ordering = false)
    @test isequal(Set(correct["nonidentifiable"]), Set(output["nonidentifiable"]))
    @test isequal(Set(correct["locally_not_globally"]), Set(output["locally_not_globally"]))
    @test isequal(Set(correct["globally"]), Set(output["globally"]))

    @test isequal(Set(output_no_weight["nonidentifiable"]), Set(output["nonidentifiable"]))
    @test isequal(Set(output_no_weight["locally_not_globally"]), Set(output["locally_not_globally"]))
    @test isequal(Set(output_no_weight["globally"]), Set(output["globally"]))
end
