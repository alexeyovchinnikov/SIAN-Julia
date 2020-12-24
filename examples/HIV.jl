using SIAN

println("Setting up the problem")

ode = @ODEmodel(
            x1'(t) = -beta * x1(t) * x4(t) - d * x1(t) + s,
            x2'(t) = beta * q1 * x1(t) * x4(t) - k1 * x2(t) - mu1 * x2(t),
            x3'(t) = beta * q2 * x1(t) * x4(t) + k1 * x2(t) - mu2 * x3(t),
            x4'(t) = -c * x4(t) + k2 * x3(t),
            y1(t) = x1(t),
            y2(t) = x4(t)
)          

identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=0, nthrds=1)


