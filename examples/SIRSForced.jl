using SIAN

println("Setting up the problem")

ode = @ODEmodel(
  s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
  i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
  r'(t) = nu * i(t) - (mu + g) * r(t),
  x1'(t) = -M * x2(t),
  x2'(t) = M * x1(t),
  y1(t) = i(t),
  y2(t) = r(t)
)          

identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=2^29 - 3, nthrds=1)


