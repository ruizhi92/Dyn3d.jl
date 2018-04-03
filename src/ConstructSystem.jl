module ConstructSystem

# export
export AddBody

# use registered packages
using DocStringExtensions

# import self-defined modules
using ..ConfigDataType

"""
This module construct the body-joint system by:
    1. AddBody
    2. AddJoint
    3. AssembleSystem
"""

# mutable struct
function AddBody(b::config_body)
    return 1
end


function AddJoint()

end

function AssembleSystem()

end






end
