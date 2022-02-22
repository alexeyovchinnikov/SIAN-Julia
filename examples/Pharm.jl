using SIAN

println("Setting up the problem")

ode = @ODEmodel(
  x1'(t) = a1 * (x2(t) - x1(t)) - (ka * n * x1(t)) / (kc * ka + kc * x3(t) + ka * x1(t)),
  x2'(t) = a1 * (x1(t) - x2(t)),
  x3'(t) = b1 * (x4(t) - x3(t)) - (kc * n * x3(t)) / (kc * ka + kc * x3(t) + ka * x1(t)),
  x4'(t) = b2 * (x3(t) - x4(t)),
  y1(t) = x1(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, nthrds = 1)

println(res)
