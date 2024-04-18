using Test
using TestSetExtensions

using SIAN
using Nemo, ModelingToolkit

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

