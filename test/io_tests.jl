@testset "Checking outputs on several examples" begin
    test_cases = []

    #---------------------------------------

    ode = @ODEmodel(
        x1'(t) = a * x1(t) - b * x1(t) * x2(t),
        x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
        y1(t) = x1(t)
    )
    correct = Dict(
        "globally" => Set([a, c, d, x1]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set([b, x2])
    )
    correct_with_known = Dict(
        "globally" => Set([x1, c, x2, d, a]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set()
    )
    known_states = [b]
    push!(
        test_cases,
        (ode, correct, correct_with_known, known_states)
    )

    #---------------------------------------

    ode = @ODEmodel(
        x1'(t) = 0,
        y1(t) = x1(t),
        y2(t) = a * x1(t) + a^2
    )
    correct = Dict(
        "globally" => Set([x1]),
        "locally_not_globally" => Set([a]),
        "nonidentifiable" => Set()
    )
    correct_with_known = Dict(
        "globally" => Set([x1]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set()
    )
    known_states = [a]
    push!(
        test_cases,
        (ode, correct, correct_with_known, known_states)
    )

    #---------------------------------------

    ode = @ODEmodel(
        S'(t) = -b * S(t) * In(t) / N(t),
        E'(t) = b * S(t) * In(t) / N(t) - nu * E(t),
        In'(t) = nu * E(t) - a * In(t),
        N'(t) = 0,
        y1(t) = In(t),
        y2(t) = N(t)
    )
    correct = Dict(
        "globally" => Set([b, In, N]),
        "locally_not_globally" => Set([nu, S, a, E]),
        "nonidentifiable" => Set()
    )
    correct_with_known = Dict(
        "globally" => Set([E, In, nu, b, N]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set()
    )
    known_states = [a, S]
    push!(
        test_cases,
        (ode, correct, correct_with_known, known_states)
    )

    #---------------------------------------

    ode = @ODEmodel(
        xA'(t) = -k1 * xA(t),
        xB'(t) = k1 * xA(t) - k2 * xB(t),
        xC'(t) = k2 * xB(t),
        eA'(t) = 0,
        eC'(t) = 0,
        y1(t) = xC(t),
        y2(t) = eA(t) * xA(t) + eB * xB(t) + eC(t) * xC(t),
        y3(t) = eA(t),
        y4(t) = eC(t)
    )
    correct = Dict(
        "globally" => Set([xC, eC, eA]),
        "locally_not_globally" => Set([xA, k1, xB, eB, k2]),
        "nonidentifiable" => Set()
    )
    correct_with_known = Dict(
        "globally" => Set([eA, xC, k2, xB, xA, eC]),
        "locally_not_globally" => Set(),
        "nonidentifiable" => Set()
    )
    known_states = [k1, eB]
    push!(
        test_cases,
        (ode, correct, correct_with_known, known_states)
    )

    #---------------------------------------

    for case in test_cases
        ode = case[1]
        result = identifiability_ode(ode, get_parameters(ode); weighted_ordering=true, known_states=[])
        result_no_weights = identifiability_ode(ode, get_parameters(ode); weighted_ordering=false, known_states=[])

        @test case[2] == result
        @test case[2] == result_no_weights

        result = identifiability_ode(ode, get_parameters(ode); weighted_ordering=true, known_states=case[4])
        result_no_weights = identifiability_ode(ode, get_parameters(ode); weighted_ordering=false, known_states=case[4])

        @test case[3] == result
        @test case[3] == result_no_weights

    end
end
