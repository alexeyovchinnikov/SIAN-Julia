using SIAN, Logging
@info "Setting up the problem"

ode = @ODEmodel(
  x₁'(t) = -1 * p₁ * x₁(t) + x₂(t) + u₀(t),
  x₂'(t) = p₃ * x₁(t) - p₄ * x₂(t) + x₃(t),
  x₃'(t) = p₆ * x₁(t) - p₇ * x₃(t),
  u₀'(t) = 1,
  y₁(t) = x₁(t),
  y₂(t) = u₀(t)
)

res = identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=0)

println(res)
