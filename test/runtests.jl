using Test
using TestSetExtensions

using SIAN
using ModelingToolkit

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

