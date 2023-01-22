using SIAN, BenchmarkTools
ode = SIAN.@ODEmodel(
    x₁'(t) = -k₁ * x₁(t) * x₂(t) + k₂ * x₃(t) + k₁₀ * x₉(t),
    x₂'(t) = -k₁ * x₁(t) * x₂(t) + k₂ * x₃(t) + k₄ * x₄(t),
    x₃'(t) = k₁ * x₁(t) * x₂(t) - (k₂ + k₃) * x₃(t),
    x₄'(t) = k₃ * x₃(t) - k₄ * x₄(t),
    x₅'(t) = k₄ * x₄(t) - k₅ * x₅(t) * x₆(t) + k₆ * x₇(t),
    x₆'(t) = -k₅ * x₅(t) * x₆(t) - k₈ * x₈(t) * x₆(t) + (k₆ + k₇) * x₇(t) + (k₉ + k₁₀) * x₉(t),
    x₇'(t) = k₅ * x₅(t) * x₆(t) - (k₆ + k₇) * x₇(t),
    x₈'(t) = k₇ * x₇(t) - k₈ * x₆(t) * x₈(t) + k₉ * x₉(t),
    x₉'(t) = k₈ * x₆(t) * x₈(t) - (k₉ + k₁₀) * x₉(t),
    y₁(t) = x₁(t),
    y₂(t) = x₂(t),
)

known_states = [x₁, x₂, x₃, x₄, x₅, x₆, x₇, x₈, x₉]
@time result = SIAN.identifiability_ode(ode; infolevel=1, local_only=false, weighted_ordering=true, p_mod=11863279, known_states=known_states)
