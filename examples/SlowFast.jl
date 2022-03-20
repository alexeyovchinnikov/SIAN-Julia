using SIAN, Logging
@info "Setting up the problem"

ode = @ODEmodel(
  xA'(t) = -k1 * xA(t),
  xB'(t) = k1 * xA(t) - k2 * xB(t),
  xC'(t) = k2 * xB(t),
  eA'(t) = 0,
  eC'(t) = 0,
  y1(t) = xC(t),
  y2(t) = eA(t) * xA(t) + eB * xB(t) + eC(t) * xC(t),
  y3(t) = eA(t),
  y4(t) = eC(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)



println(res)

# The following computation shows that one can identify more from two experiments
# SIAN. generate_replica(ode, 2) generates a system consiting of two copies of sigma
# with the same parameters in ODEs but different inputs and initial conditions

r = 2

res = identifiability_ode(SIAN.generate_replica(ode, r), get_parameters(ode; initial_conditions = false); p = 0.99, p_mod = 0, nthrds = 1)

println(res)
