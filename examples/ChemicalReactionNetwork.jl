using Pkg; Pkg.activate("../IdentifiabilityODE"); using IdentifiabilityODE

println("Setting up the problem")

ode = @ODEmodel(
  x1'(t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k4 * x6(t),
  x2'(t) = k1 * x1(t) * x2(t) + k2 * x4(t) + k3 * x4(t),
  x3'(t) = k3 * x4(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
  x4'(t) = k1 * x1(t) * x2(t) - k2 * x4(t) - k3 * x4(t),
  x5'(t) = k4 * x6(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
  x6'(t) = -k4 * x6(t) - k5 * x6(t) + k6 * x3(t) * x5(t),
  y1(t) = x3(t),
  y2(t) = x2(t)
)

res = identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=2^29 - 3, nthrds=3)

println(res)
