using Test
using TestSetExtensions

using Pkg; Pkg.activate("../IdentifiabilityODE"); using IdentifiabilityODE

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

