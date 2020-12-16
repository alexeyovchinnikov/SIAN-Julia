using Test
using TestSetExtensions

using SIAN

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

