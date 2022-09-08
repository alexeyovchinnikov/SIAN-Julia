# BIOMD0000000003
using SIAN, BenchmarkTools

ode = @ODEmodel(
    x1'(t) = (1 * k6 * k7 + (-1) * x1(t) * k6 * k8 + (-1) * x1(t) * k6 * k9 * x3(t) * 1 / (x1(t) + k10)) / k6,
    x2'(t) = (1 * k6 * (1 + (-1) * x2(t)) * x1(t) * k3 * 1 / (x1(t) + k5) * 1 / (k11 + (-1) * x2(t) + 1) + (-1) * k6 * x2(t) * k12 * 1 / (k13 + x2(t))) / k6,
    x3'(t) = (1 * k6 * x2(t) * k4 * (1 + (-1) * x3(t)) * 1 / (k14 + (-1) * x3(t) + 1) + (-1) * k6 * k16 * x3(t) * 1 / (k15 + x3(t))) / k6,
    y1(t) = x1(t)
)

known_states = [x1]

@time result = identifiability_ode(ode; infolevel=0, local_only=false, weighted_ordering=true, known_states=[])
@time result = identifiability_ode(ode; infolevel=0, local_only=false, weighted_ordering=true, known_states=known_states)