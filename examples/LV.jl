using Pkg; Pkg.activate("../IdentifiabilityODE"); using IdentifiabilityODE

println("Setting up the problem")

ode = @ODEmodel(
    x1'(t) = a * x1(t) - b * x1(t) * x2(t),
    x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
    y1(t) = x1(t)
)

res = identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=0, nthrds=1)

println(res)


