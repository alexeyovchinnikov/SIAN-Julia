@testset "Testing weight generation" begin
    # test actual array of weights
    @info "Lotka-Volterra"
    ode = @ODEmodel(
        x1'(t) = a * x1(t) - b * x1(t) * x2(t),
        x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
        y1(t) = x1(t)
    )

    param2str = Dict(string(p) => p for p in vcat(ode.parameters, ode.x_vars, ode.y_vars))
    maple_weights = Dict("c" => 3, "x2" => 2, "x1" => 1, "z_aux" => 2, "d" => 3)
    result = SIAN.identifiability_ode(ode, SIAN.get_parameters(ode); local_only = true)

    weights = SIAN.get_weights(ode, result["non_identifiable"])
    @test all(maple_weights[string(k)] == weights[k] for k in keys(weights))

    @info "Cholera"
    ode = @ODEmodel(
        s'(t) = mu - bi * s(t) * i(t) - bw * s(t) * w(t) - mu * s(t) + al * r(t),
        i'(t) = bw * s(t) * w(t) + bi * s(t) * i(t) - g * i(t) - mu * i(t),
        w'(t) = dz * (i(t) - w(t)),
        r'(t) = g * i(t) - mu * r(t) - al * r(t),
        y1(t) = k * i(t),
        y2(t) = i(t) + r(t) + s(t)
    )

    param2str = Dict(string(p) => p for p in vcat(ode.parameters, ode.x_vars, ode.y_vars))
    maple_weights = Dict("dz" => 3, "z_aux" => 2, "r" => 1, "w" => 2, "i" => 1, "s" => 1, "al" => 3)

    res = identifiability_ode(ode, get_parameters(ode); local_only = true)
    weights = SIAN.get_weights(ode, result["non_identifiable"])
    @test all(maple_weights[string(k)] == weights[k] for k in keys(weights))

    @info "SIR_R₀"
    ode = @ODEmodel(
        S'(t) = -b * In(t) * S(t),
        In'(t) = b * In(t) * S(t) - g * In(t),
        R'(t) = g * In(t),
        aux'(t) = 0,
        y1(t) = In(t),
        y2(t) = b // g + aux(t)
    )
    param2str = Dict(string(p) => p for p in vcat(ode.parameters, ode.x_vars, ode.y_vars))
    maple_weights = Dict("z_aux" => 1, "In" => 1, "S" => 2)

    res = identifiability_ode(ode, get_parameters(ode); local_only = true)
    weights = SIAN.get_weights(ode, result["non_identifiable"])
    @test all(maple_weights[string(k)] == weights[k] for k in keys(weights))

    @info "Lipolysis"

    ode = @ODEmodel(
        x1'(t) = -x1(t) * x5(t) / (k2 + x1(t)),
        x2'(t) = 2 * x1(t) * x5(t) / ((k2 + x1(t)) * 3) - k4 * x2(t),
        x3'(t) = k4 * (x2(t)) / 2 - k4 * x3(t),
        x4'(t) = x1(t) * x5(t) / (3 * (k2 + x1(t))) + k4 * (x2(t)) / 2 + k4 * x3(t),
        x5'(t) = -k3 * x5(t),
        y1(t) = x1(t),
        y2(t) = x2(t) + x3(t),
        y3(t) = x4(t)
    )
    param2str = Dict(string(p) => p for p in vcat(ode.parameters, ode.x_vars, ode.y_vars))
    maple_weights = Dict("k3" => 3, "x2" => 1, "x1" => 1, "z_aux" = 2, "x3" => 1, "x5" => 2)

    res = identifiability_ode(ode, get_parameters(ode); local_only = true)
    weights = SIAN.get_weights(ode, result["non_identifiable"])
    @test all(maple_weights[string(k)] == weights[k] for k in keys(weights))
end
