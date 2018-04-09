module ConstructSystem

# export
export SingleBody, SingleJoint, AddBody, AddJoint, AssembleSystem

# use registered packages
using DocStringExtensions
import Base: show

# import self-defined modules
using ..ConfigDataType

"""
This module construct the body-joint system by:
    1. AddBody
    2. AddJoint
    3. AssembleSystem
"""

#-------------------------------------------------------------------------------
mutable struct SingleBody
    # hierarchy info
    id::Int
    pid::Int
    chid::Vector{Int}
    nchild::Int
    # verts
    nverts::Int
    verts::Array{Float64,2}
    verts_i::Array{Float64,2}
    # coord in body frame
    x_c::Vector{Float64}
    x_0::Vector{Float64}
    # mass and inertia
    mass::Float64
    inertia_c::Array{Float64,2}
    inertia_0::Array{Float64,2}
    # transform matrix
    Xb_to_c::Array{Float64,2}
    Xb_to_i::Array{Float64,2}
    Xp_to_b::Array{Float64,2}
    # motion vectors
    q::Vector{Float64}
    v::Vector{Float64}
    c::Vector{Float64}
    # articulated body info
    pA::Vector{Float64}
    Ib_A::Array{Float64,2}
end

# outer constructor
SingleBody() = SingleBody(
    0,0,Vector{Int}(0),0,
    0,Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    Vector{Float64}(0),Vector{Float64}(0),
    0,Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    Array{Float64,2}(0,0),Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    Vector{Float64}(0),Vector{Float64}(0),Vector{Float64}(0),
    Vector{Float64}(0),Array{Float64,2}(0,0)
)

function show(io::IO, m::SingleBody)
    println(io, "body_id = $(m.id)", ", parent_id = $(m.pid)",
            ", nchild = $(m.nchild)", ", child_id = $(m.chid)")
    println(io, "nverts = $(m.nverts)", ", verts = $(m.verts)")
    println(io, "verts_i = $(m.verts_i)")
    println(io, "x_c = $(m.x_c)", ", x_0 = $(m.x_0)")
    println(io, "mass = $(m.mass)")
    println(io, "inertia_c = $(m.inertia_c)")
    println(io, "inertia_0 = $(m.inertia_0)")
    println(io, "Xb_to_c = $(m.Xb_to_c)")
    println(io, "Xb_to_i = $(m.Xb_to_i)")
    println(io, "Xp_to_b = $(m.Xp_to_b)")
    println(io, "q = $(m.q)", ", v = $(m.v)", ", c = $(m.c)")
end

#-------------------------------------------------------------------------------
mutable struct SingleJoint


end

#-------------------------------------------------------------------------------
# mutable struct
function AddBody(b::ConfigBody)
    return 1
end


function AddJoint()

end

function AssembleSystem()

end






end
