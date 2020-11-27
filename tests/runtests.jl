using Test
using TestSetExtensions

include("../IdentifiabilityODE.jl")

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

