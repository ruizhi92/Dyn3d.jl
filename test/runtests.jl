using Base.Test
using TestSetExtensions

include(Pkg.dir("Dyn3d")*"/src/Dyn3d.jl")
using Dyn3d
include(Pkg.dir("Whirl")*"/src/Whirl.jl")
using Whirl

@test isempty(detect_ambiguities(Dyn3d))

@testset ExtendedTestSet "All tests" begin
    @includetests ARGS
end

if isempty(ARGS)
    println("Test done.")
end
