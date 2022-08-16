using SIAN, Logging
@info "Setting up the problem"

ode = @ODEmodel(
    x₁'(t) = -β * x₁(t) * x₄(t) - d * x₁(t) + s,
    x₂'(t) = β * q₁ * x₁(t) * x₄(t) - k₁ * x₂(t) - μ₁ * x₂(t),
    x₃'(t) = β * q₂ * x₁(t) * x₄(t) + k₁ * x₂(t) - μ₂ * x₃(t),
    x₄'(t) = -c * x₄(t) + k₂ * x₃(t),
    y₁(t) = x₁(t),
    y₂(t) = x₄(t)
)

res = identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=0)

println(res)
