println("Setting up the problem")

include("IdentifiabilityODE.jl")

using BenchmarkTools

prob=0.99
R, (x1,x2,y1,u,v,a,b,c,d) = Nemo.PolynomialRing(Nemo.QQ, ["x1","x2","y1","u","v","a","b","c","d"], ordering=:degrevlex)

sigma_x = [
            [x1, a*x1 - b*x1*x2],
            [x2, -c*x2 + d*x1*x2],
          ]
sigma_y = [
            [y1, x1+u+v]
          ]

r=3

sigma_x_replica, sigma_y_replica, u_vars_replica = GenerateReplica(sigma_x, sigma_y,  [u,v], r)

for i in 1:r*length(sigma_x)
   println(sigma_x_replica[i])
end

for i in 1:r*length(sigma_y)             
   println(sigma_y_replica[i])
end

println(u_vars_replica)


IdentifiabilityODE(sigma_x, sigma_y,  [u,v], GetParameters(sigma_x, sigma_y, [u,v]), prob)


