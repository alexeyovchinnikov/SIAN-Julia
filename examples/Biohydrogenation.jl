using SIAN, Logging

@info "Setting up the problem"

ode = @ODEmodel(
  x₄'(t) = -k₅ * x₄(t) // (k₆ + x₄(t)),
  x₅'(t) = k₅ * x₄(t) // (k₆ + x₄(t)) - k₇ * x₅(t) / (k₈ + x₅(t) + x₆(t)),
  x₆'(t) = k₇ * x₅(t) // (k₈ + x₅(t) + x₆(t)) - k₉ * x₆(t) * (k₁₀ - x₆(t)) // k₁₀,
  x₇'(t) = k₉ * x₆(t) * (k₁₀ - x₆(t)) // k₁₀,
  y₁(t) = x₄(t),
  y₂(t) = x₅(t)
)

res = identifiability_ode(ode, get_parameters(ode); p=0.999, p_mod=0)

println(res)
