module UpdateSystem

export UpdateVelocity, UpdatePosition

# use registered packages
using DocStringExtensions

# abstract type Velocity end
# abstract type Position end

#-------------------------------------------------------------------------------
function UpdatePosition(kind::String, qJ::Vector{T}) where T <: AbstractFloat
    if kind == "revolute" || kind == "prismatic" || kind == "cylindrical" ||
       kind == "planar" || kind == "extended_revolute"
       X = TransMatrix(qJ)
   elseif kind == "helical"
       # set fixed screw parameter h
       h = 0.1
       r = h*qJ[4:6]
       θ = qJ[1:3]
       X = TransMatrix([θ; r])
   elseif kind == "spherical"
       error("Under construction")
   elseif kind == "free"
       error("Under construction")
   end
   return X
end

#-------------------------------------------------------------------------------
function UpdateVelocity()
"""
"""
end



end
