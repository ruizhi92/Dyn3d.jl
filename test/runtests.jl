using Base.Test
using TestSetExtensions

include("../src/Dyn3d.jl")
using Dyn3d

@test isempty(detect_ambiguities(Dyn3d))

@testset ExtendedTestSet "All tests" begin
    @includetests ARGS
end

if isempty(ARGS)
    println("Test done.")
end
