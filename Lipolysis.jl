include("IdentifiabilityODE.jl")

println("Setting up the problem")

ode = @ODEmodel(
  x1'(t) = -x1(t) * x5(t) / (k2 + x1(t)),
  x2'(t) = 2 * x1(t) * x5(t) / ((k2 + x1(t)) * 3) - k4 * x2(t),
  x3'(t) = k4*(x2(t))/2 - k4*x3(t),
  x4'(t) = x1(t) * x5(t) / (3 * (k2 + x1(t))) + k4 * (x2(t))/2 + k4 * x3(t),
  x5'(t) = -k3 * x5(t),
  y1(t) = x1(t),
  y2(t) = x2(t) + x3(t),
  y3(t) = x4(t)
)          

identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 64)


