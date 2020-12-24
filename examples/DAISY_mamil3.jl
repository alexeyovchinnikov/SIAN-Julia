using SIAN

println("Setting up the problem")

ode = @ODEmodel(
  x1'(t) = -(a21 + a31 + a01) * x1(t) + a12 * x2(t) + a13 * x3(t) + u(t),
  x2'(t) = a21 * x1(t) - a12 * x2(t),
  x3'(t) = a31 * x1(t) - a13 * x3(t),
  y(t) = x1(t)
)          

identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=0, nthrds=1)


