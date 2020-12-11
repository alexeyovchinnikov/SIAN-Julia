using Test
using TestSetExtensions

include("../src/IdentifiabilityODE.jl")

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

