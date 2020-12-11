include("../src/IdentifiabilityODE.jl")

println("Setting up the problem")

ode = @ODEmodel(
  x'(t) = lm - d * x(t) - beta * x(t) * v(t),
  y'(t) = beta * x(t) * v(t) - a * y(t),
  v'(t) = k * y(t) - u * v(t),
  w'(t) = c * z(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
  z'(t) = c * q * y(t) * w(t) - h * z(t),
  y1(t) = w(t),
  y2(t) = z(t)
)          

identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, nthrds = 64)


