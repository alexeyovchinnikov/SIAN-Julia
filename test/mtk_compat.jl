@testset "Check identifiability of `ODESystem` object" begin
    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y1(t)
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1]
    measured_quantities = [y1 ~ x0]
    de = ODESystem(eqs, t, name=:Test)
    correct = Dict{Any,Any}(
        "nonidentifiable" => [a12, a21, substitute(x1, t => 0), a01],
        "locally_not_globally" => Num[],
        "globally" => [substitute(x0, t => 0)]
    )
    output = identifiability_ode(de; measured_quantities=measured_quantities)

    @test isequal(Set(correct["nonidentifiable"]), Set(output["nonidentifiable"]))
    @test isequal(Set(correct["locally_not_globally"]), Set(output["locally_not_globally"]))
    @test isequal(Set(correct["globally"]), Set(output["globally"]))

    # ---------------------------------------

    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y1(t)
    D = Differential(t)
    using SpecialFunctions

    eqs = [
        D(x0) ~ -(a01 + a21) * SpecialFunctions.erfc(x0) + a12 * x1,
        D(x1) ~ a21 * x0 - a12 * x1
    ]
    measured_quantities = [y1 ~ x0]

    de = ODESystem(eqs, t, name=:Test)

    @test_throws ArgumentError identifiability_ode(de; measured_quantities=measured_quantities)
    # ---------------------------------------

    test_cases = []
    # ---------------------------------------
    @parameters a, b, c, d
    @variables t x1(t) x2(t) y1(t)

    eqs = [
        D(x1) ~ a * x1 - b * x1 * x2,
        D(x2) ~ -c * x2 + d * x1 * x2,
    ]
    measured_quantities = [y1 ~ x1]
    known_states = [b]
    @named ode = ODESystem(eqs, t)

    correct = Dict(
        "globally" => Set([a, c, d, substitute(x1, t => 0)]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set([b, substitute(x2, t => 0)])
    )
    correct_with_known = Dict(
        "globally" => Set([substitute(x1, t => 0), c, substitute(x2, t => 0), d, a]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set()
    )

    push!(
        test_cases,
        (ode, correct, correct_with_known, known_states, measured_quantities)
    )

    #---------------------------------------

    @parameters a
    @variables t x1(t) y1(t) y2(t)
    eqs = [
        D(x1) ~ a * x1,
    ]
    measured_quantities = [y1 ~ x1, y2 ~ a * x1 + a^2]
    known_states = [a]
    @named ode = ODESystem(eqs, t)
    correct = Dict(
        "globally" => Set([substitute(x1, t => 0)]),
        "locally_not_globally" => Set([a]),
        "nonidentifiable" => Set()
    )
    correct_with_known = Dict(
        "globally" => Set([substitute(x1, t => 0)]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set()
    )

    push!(
        test_cases,
        (ode, correct, correct_with_known, known_states, measured_quantities)
    )

    #---------------------------------------
    @parameters b nu a
    @variables t S(t) E(t) In(t) N(t) y1(t) y2(t)
    eqs = [
        D(S) ~ -b * S * In / N,
        D(E) ~ b * S * In / N - nu * E,
        D(In) ~ nu * E - a * In,
        D(N) ~ 0,
    ]
    measured_quantities = [y1 ~ In, y2 ~ N]
    known_states = [a, S]
    @named ode = ODESystem(eqs, t)

    correct = Dict(
        "globally" => Set([b, substitute(In, t => 0), substitute(N, t => 0)]),
        "locally_not_globally" => Set([nu, substitute(S, t => 0), a, substitute(E, t => 0)]),
        "nonidentifiable" => Set()
    )
    correct_with_known = Dict(
        "globally" => Set([substitute(E, t => 0), substitute(In, t => 0), nu, b, substitute(N, t => 0)]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set()
    )

    push!(
        test_cases,
        (ode, correct, correct_with_known, known_states, measured_quantities)
    )

    #---------------------------------------
    @parameters k1 k2 eB
    @variables t xA(t) xB(t) xC(t) eA(t) eC(t) y1(t) y2(t) y3(t) y4(t)

    eqs = [
        D(xA) ~ -k1 * xA,
        D(xB) ~ k1 * xA - k2 * xB,
        D(xC) ~ k2 * xB,
        D(eA) ~ 0,
        D(eC) ~ 0,
    ]
    measured_quantities = [y1 ~ xC, y2 ~ eA * xA + eB * xB + eC * xC, y3 ~ eA, y4 ~ eC]
    known_states = [k1, eB]
    @named ode = ODESystem(eqs, t)

    correct = Dict(
        "globally" => Set([substitute(xC, t => 0), substitute(eC, t => 0), substitute(eA, t => 0)]),
        "locally_not_globally" => Set([substitute(xA, t => 0), k1, substitute(xB, t => 0), eB, k2]),
        "nonidentifiable" => Set()
    )
    correct_with_known = Dict(
        "globally" => Set([substitute(eA, t => 0), substitute(xC, t => 0), k2, substitute(xB, t => 0), substitute(xA, t => 0), substitute(eC, t => 0)]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set()
    )

    push!(
        test_cases,
        (ode, correct, correct_with_known, known_states, measured_quantities)
    )

    #---------------------------------------

    for case in test_cases
        ode = case[1]
        result = identifiability_ode(ode; weighted_ordering=true, known_states=[])
        result_no_weights = identifiability_ode(ode; measured_quantities=case[end], weighted_ordering=false, known_states=[])

        @test case[2] == result
        @test case[2] == result_no_weights

        result = identifiability_ode(ode; measured_quantities=case[end], weighted_ordering=true, known_states=case[4])
        result_no_weights = identifiability_ode(ode; measured_quantities=case[end], weighted_ordering=false, known_states=case[4])

        @test case[3] == result
        @test case[3] == result_no_weights

    end
end
