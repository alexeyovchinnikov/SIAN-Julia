println("Setting up the problem")

include("IdentifiabilityODE.jl")

prob=0.99
R, (x1,x2,y1,u,a,b,c,d) = Oscar.PolynomialRing(Oscar.QQ, ["x1","x2","y1","u","a","b","c","d"], ordering=:degrevlex)

sigma_x = [
            [x1, a*x1 - b*x1*x2],
            [x2, -c*x2 + d*x1*x2],
          ]
sigma_y = [
            [y1, x1+u]
          ]

IdentifiabilityODE(sigma_x, sigma_y,  [u], [a,b,c,d], [x1], prob)


