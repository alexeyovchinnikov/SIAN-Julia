include("IdentifiabilityODE.jl")

println("Setting up the problem")

R, (x1, x2, x3, x4, y1, a1, b1, b2, ka, n, kc) = Nemo.PolynomialRing(Nemo.QQ, ["x1", "x2", "x3", "x4", "y1", "a1", "b1", "b2", "ka", "n", "kc"])

sigma_x = [
            [x1, a1 * (x2 - x1) - (ka * n * x1) // (kc * ka + kc * x3 + ka * x1)],
            [x2, a1 * (x1 - x2)],
            [x3, b1 * (x4 - x3) - (kc * n * x3) // (kc * ka + kc * x3 + ka * x1)],
            [x4, b2 * (x3 - x4)]
          ]
sigma_y = [
            [y1, x1]
          ]

identifiability_ode(sigma_x, sigma_y,  [], get_parameters(sigma_x, sigma_y, []); p = 0.99, p_mod = 2^29 - 3, nthrds = 64)


