# Taken from 
# Vajda S., Rabitz H.
# Identifiability and Distinguishability of First-Order Reaction Systems, p. 701
# We added an extra output x_C
read "IdentifiabilityODE.mpl":

sigma := [
  diff(xA(t), t) = -k1 * xA(t),
  diff(xB(t), t) = k1 * xA(t) - k2 * xB(t),
  diff(xC(t), t) = k2 * xB(t),
  diff(eA(t), t) = 0,
  diff(eC(t), t) = 0,
  y1(t) = xC(t),
  y2(t) = eA(t) * xA(t) + eB * xB(t) + eC(t) * xC(t),
  y3(t) = eA(t),
  y4(t) = eC(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma));

# The following computation shows that one can identify more from two experiments
# GenerateReplica(sigma, 2) generates a system consiting of two copies of sigma
# with the same parameters in ODEs but different inputs and initial conditions

IdentifiabilityODE(GenerateReplica(sigma, 2), GetParameters(sigma, initial_conditions=false));
