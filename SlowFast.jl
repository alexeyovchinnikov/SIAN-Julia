println("Setting up the problem")

include("IdentifiabilityODE.jl")

prob=0.99

R, (xA, xB, xC, eA, eB, eC, y1, y2, y3, y4, k1, k2) = Nemo.PolynomialRing(Nemo.QQ, ["xA", "xB", "xC", "eA", "eB", "eC", "y1", "y2", "y3", "y4", "k1", "k2"])

sigma_x = [
            [xA, -k1 * xA],
            [xB, k1 * xA - k2 * xB],
            [xC, k2 * xB],
            [eA, R(0)],
            [eC, R(0)]
          ]
sigma_y = [
            [y1, xC],
            [y2, eA * xA + eB * xB + eC * xC],
            [y3, eA],
            [y4, eC]
          ]

identifiability_ode(sigma_x, sigma_y, [], get_parameters(sigma_x, sigma_y, []), prob);

# The following computation shows that one can identify more from two experiments
# GenerateReplica(sigma, 2) generates a system consiting of two copies of sigma
# with the same parameters in ODEs but different inputs and initial conditions

r=2

sigma_x_r, sigma_y_r, u_vars_r = generate_replica(sigma_x, sigma_y, [], r)

identifiability_ode(sigma_x_r, sigma_y_r, u_vars_r, get_parameters(sigma_x_r, sigma_y_r, u_vars_r; initial_conditions = false), prob)


