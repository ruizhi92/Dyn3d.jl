__precompile__()

module Dyn3d

using Reexport
using DocStringExtensions

include("config_files/ConfigDataType.jl")
@reexport using .ConfigDataType

include("Utils.jl")
@reexport using .Utils

include("SpatialAlgebra.jl")
@reexport using .SpatialAlgebra

include("ConstructSystem.jl")
@reexport using .ConstructSystem

include("UpdateSystem.jl")
@reexport using .UpdateSystem

include("TimeMarching.jl")
@reexport using .TimeMarching

#=
21-element Array{Coverage.MallocInfo,1}:
 Coverage.MallocInfo(8658000, "../src/TimeMarching.jl.mem", 230)
 Coverage.MallocInfo(8658000, "../src/TimeMarching.jl.mem", 279)
 Coverage.MallocInfo(10158720, "../src/TimeMarching.jl.mem", 200)
 Coverage.MallocInfo(10620480, "../src/TimeMarching.jl.mem", 197)
 Coverage.MallocInfo(11544000, "../src/TimeMarching.jl.mem", 226)
 Coverage.MallocInfo(11544000, "../src/TimeMarching.jl.mem", 250)
 Coverage.MallocInfo(11544000, "../src/TimeMarching.jl.mem", 275)
 Coverage.MallocInfo(12931072, "../src/UpdateSystem.jl.mem", 84)
 Coverage.MallocInfo(13391040, "../src/TimeMarching.jl.mem", 202)
 Coverage.MallocInfo(14778368, "../src/UpdateSystem.jl.mem", 81)
 Coverage.MallocInfo(14780416, "../src/SpatialAlgebra.jl.mem", 36)
 Coverage.MallocInfo(17316000, "../src/TimeMarching.jl.mem", 254)
 Coverage.MallocInfo(18011136, "../src/UpdateSystem.jl.mem", 105)
 Coverage.MallocInfo(27713280, "../src/SpatialAlgebra.jl.mem", 22)
 Coverage.MallocInfo(32096768, "../src/UpdateSystem.jl.mem", 70)
 Coverage.MallocInfo(47099520, "../src/TimeMarching.jl.mem", 205)
 Coverage.MallocInfo(60028800, "../src/TimeMarching.jl.mem", 303)
 Coverage.MallocInfo(96268608, "../src/UpdateSystem.jl.mem", 95)
 Coverage.MallocInfo(125367840, "../src/TimeMarching.jl.mem", 144)
 Coverage.MallocInfo(131176192, "../src/SpatialAlgebra.jl.mem", 39)
 Coverage.MallocInfo(146839680, "../src/TimeMarching.jl.mem", 134)
=#

end
