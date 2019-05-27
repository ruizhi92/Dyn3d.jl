using Pkg
Pkg.activate("..")

using Compat
using Compat.Test

using Dyn3d

include("config_body.jl")
include("ConstructSystem.jl")
include("Timemarching.jl")
include("SpatialAlgebra.jl")
include("FluidInteraction.jl")
include("Utils.jl")
