include("IdentifiabilityODE.jl")

println("Setting up the problem")

R, (x1, x2, y1, u, v, a, b, c, d) = Nemo.PolynomialRing(Nemo.QQ, ["x1", "x2", "y1", "u", "v", "a", "b", "c", "d"])

sigma_x = [
            [x1, a * x1 - b * x1 * x2],
            [x2, -c * x2 + d * x1 * x2],
          ]
sigma_y = [
            [y1, x1 + u + v]
          ]

identifiability_ode(sigma_x, sigma_y,  [u, v], get_parameters(sigma_x, sigma_y, [u,v]); p = 0.99, p_mod = 2^29 - 3, nthrds = 64)


