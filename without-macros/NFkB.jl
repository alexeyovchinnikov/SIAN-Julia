include("IdentifiabilityODE.jl")

println("Setting up the problem")

R, (x1, x2, x3, x4, x5, x6, x7, x9, x10, x11, x12, x13, x14, x15, y1, y2, y3, y4, y5, x8, y6, u, c5, c_3a, c_4a, e_2a, i1, i_1a, k1, k2, k3, k_deg, k_prod, t1, t2) = Nemo.PolynomialRing(Nemo.QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "y1", "y2", "y3", "y4", "y5", "x8", "y6", "u", "c5", "c_3a", "c_4a", "e_2a", "i1", "i_1a", "k1", "k2", "k3", "k_deg", "k_prod", "t1", "t2"])

a1 = R(1 // 2)
a2 = R(1 // 5)
a3 = R(1)
c_1a = R(5 // 10^(7))
c_2a = R(0)
c_5a = R(1 // 10^(4))
c_6a = R(2 // 10^(5))
c1 = R(5 // 10^(7))
c2 = R(0)
c3 = R(4 // 10^(4))
c4 = R(1 // 2)
kv = R(5)
e_1a = R(5 // 10^(4))
c_1c = R(5 // 10^(7))
c_2c = R(0)
c_3c = R(4 // 10^(4))

sigma_x = [
  [x1, k_prod - k_deg * x1 - k1 * x1 * u],
  [x2, -k3 * x2 - k_deg * x2 - a2 * x2 * x10 + t1 * x4 - a3 * x2 * x13 + t2 * x5 + (k1 * x1 - k2 * x2 * x8) * u],
  [x3,  k3 * x2 - k_deg * x3 + k2 * x2 * x8 * u],
  [x4, a2 * x2 * x10 - t1 * x4],
  [x5, a3 * x2 * x13 - t2 * x5],
  [x6, c_6a * x13 - a1 * x6 * x10 + t2 * x5 - i1 * x6],
  [x7, i1 * kv * x6 - a1 * x11 * x7],
  [x8, c4 * x9 - c5 * x8],
  [x9, c2 + c1 * x7 - c3 * x9],
  [x10, -a2 * x2 * x10 - a1 * x10 * x6 + c_4a * x12 - c_5a * x10 - i_1a * x10 + e_1a * x11],
  [x11, -a1 * x11 * x7 + i_1a * kv * x10 - e_1a * kv * x11],
  [x12,  c_2a + c_1a * x7 - c_3a * x12],
  [x13,   a1 * x10 * x6 - c_6a * x13 - a3 * x2 * x13 + e_2a * x14],
  [x14,  a1 * x11 * x7 - e_2a * kv * x14],
  [x15,  c_2c + c_1c * x7 - c_3c * x15]
]  

sigma_y= [
  [y1, x2],
  [y2, x10 + x13],
  [y3, x9],
  [y4, x1 + x2 + x3],
  [y6, x12],
  [y5, x7]
]

identifiability_ode(sigma_x, sigma_y,  [u], get_parameters(sigma_x, sigma_y, [u]); p = 0.99, p_mod = 2^29 - 3, nthrds = 64)
