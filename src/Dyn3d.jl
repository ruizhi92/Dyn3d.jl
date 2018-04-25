module Dyn3d
#=
    Include this module for any examples in notebook.
=#

include("Config_files/ConfigDataType.jl")
using .ConfigDataType

include("SpatialAlgebra.jl")
using .SpatialAlgebra

include("ConstructSystem.jl")
using .ConstructSystem

include("UpdateSystem.jl")
using .UpdateSystem

include("TimeMarching.jl")
using .TimeMarching

include("Utils.jl")
using .Utils



end
