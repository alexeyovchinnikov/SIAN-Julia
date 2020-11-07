include("../IdentifiabilityODE.jl")

println("Setting up the problem")

ode = @ODEmodel(
  x1'(t) = -(k3 + k7) * x1(t) + k4 * x2(t),
  x2'(t) = k3 * x1(t) - (k4 + a(t) * k5 + b(t) * d(t) * k5) * x2(t) + k6 * x3(t) + k6 * x4(t) + k5 * x2(t) * x3(t) + k5 * x2(t) * x4(t),
  x3'(t) = a(t) * k5 * x2(t) - k6 * x3(t) - k5 * x2(t) * x3(t),
  x4'(t) = b(t) * d(t) * k5 * x2(t) - k6 * x4(t) - k5 * x2(t) * x4(t),
  x5'(t) = k7 * x1(t),
  a'(t) = 0,
  b'(t) = 0,
  d'(t) = 0,
  y1(t) = x5(t),
  y2(t) = a(t),
  y3(t) = b(t),
  y4(t) = d(t)
)          

identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, nthrds = 64)


