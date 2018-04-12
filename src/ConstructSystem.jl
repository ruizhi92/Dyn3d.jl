module ConstructSystem

# export
export SingleBody, SingleJoint, AddBody, AddJoint, AssembleSystem

# use registered packages
using DocStringExtensions
import Base: show

# import self-defined modules
using ..ConfigDataType
using ..SpatialAlgebra

"""
This module construct the body-joint system by:
    1. AddBody
    2. AddJoint
    3. AssembleSystem
"""

#-------------------------------------------------------------------------------
mutable struct SingleBody
    # hierarchy info
    bid::Int
    pid::Int
    chid::Vector{Int}
    nchild::Int
    # verts
    nverts::Int
    verts::Array{Float64,2}
    verts_i::Array{Float64,2}
    # coord in body frame and inertial frame
    x_c::Vector{Float64}
    x_i::Vector{Float64}
    # mass and inertia
    mass::Float64
    inertia_c::Array{Float64,2}
    inertia_i::Array{Float64,2}
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

function show(io::IO, ::MIME"text/plain", m::SingleBody)
    println(io, "body_id = $(m.bid)", ", parent_id = $(m.pid)",
            ", nchild = $(m.nchild)", ", child_id = $(m.chid)")
    println(io, "nverts = $(m.nverts)", ", verts = $(m.verts)")
    println(io, "verts_i = $(m.verts_i)")
    println(io, "x_c = $(m.x_c)", ", x_i = $(m.x_i)")
    println(io, "mass = $(m.mass)")
    println(io, "inertia_c = $(m.inertia_c)")
    println(io, "inertia_i = $(m.inertia_i)")
    println(io, "Xb_to_c = $(m.Xb_to_c)")
    println(io, "Xb_to_i = $(m.Xb_to_i)")
    println(io, "Xp_to_b = $(m.Xp_to_b)")
    println(io, "q = $(m.q)")
    println(io, "v = $(m.v)")
    println(io, "c = $(m.c)")
end

#-------------------------------------------------------------------------------
mutable struct SingleJoint
    # hierarchy info
    jid::Int
    joint_type::String
    body1::Int
    shape1::Vector{Float64}
    shape2::Vector{Float64}
    # dof info
    nudof::Int
    ncdof::Int
    np::Int
    na::Int
    udof::Vector{Float64} # newline in constructor
    cdof::Vector{Float64}
    udof_p::Vector{Float64} # newline in constructor
    udof_a::Vector{Float64}
    i_udof_p::Vector{Float64} # newline in constructor
    i_udof_a::Vector{Float64}
    udof_map::Vector{Float64}
    # dof info modified by HERK
    nudof_HERK::Int
    ncdof_HERK::Int
    udof_HERK::Vector{Float64}
    cdof_HERK::Vector{Float64}
    # joint basis matrix
    S::Array{Float64,2}
    T::Array{Float64,2}
    T_HERK::Array{Float64,2}
    # joint dof
    joint_dof::Vector{Dof}
    # motion vectors
    qJ::Vector{Float64}
    vJ::Vector{Float64}
    cJ::Vector{Float64}
    # transform matrix
    Xj::Array{Float64,2}
    Xp_to_j::Array{Float64,2}
    Xj_to_ch::Array{Float64,2}
end

# outer constructor
SingleJoint() = SingleJoint(
    0," ",0,Vector{Float64}(0),Vector{Float64}(0),
    0,0,0,0,
    Vector{Float64}(0),Vector{Float64}(0),
    Vector{Float64}(0),Vector{Float64}(0),
    Vector{Float64}(0),Vector{Float64}(0),Vector{Float64}(0),
    0,0,Vector{Float64}(0),Vector{Float64}(0),
    Array{Float64,2}(0,0),Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    [Dof()], # dof
    Vector{Float64}(0),Vector{Float64}(0),Vector{Float64}(0),
    Array{Float64,2}(0,0),Array{Float64,2}(0,0),Array{Float64,2}(0,0)
)

function show(io::IO, m::SingleJoint)
    println(io, "joint_id = $(m.jid)", ", joint_type = $(m.joint_type)",
            ", pid = $(m.body1)")
    println(io, "shape1 = $(m.shape1)", ", shape2 = $(m.shape2)")
    println(io, "nudof = $(m.nudof)", ", ncdof = $(m.ncdof)",
            ", np = $(m.np)", ", na = $(m.na)")
    println(io, "udof = $(m.udof)", ", cdof = $(m.cdof)")
    println(io, "udof_p = $(m.udof_p)", ", udof_a = $(m.udof_a)")
    println(io, "i_udof_p = $(m.i_udof_p)", ", i_udof_a = $(m.i_udof_a)")
    println(io, "udof_map = $(m.udof_map)")
    println(io, "nudof_HERK = $(m.nudof_HERK)",", ncdof_HERK = $(m.ncdof_HERK)")
    println(io, "udof_HERK = $(m.udof_HERK)", ", cdof_HERK = $(m.cdof_HERK)")
    println(io, "S = $(m.S)")
    println(io, "T = $(m.T)")
    println(io, "T_HERK = $(m.T_HERK)")
    println(io, "joint_dof = $(m.joint_dof)")
    println(io, "Xj = $(m.Xj)")
    println(io, "Xp_to_j = $(m.Xp_to_j)")
    println(io, "Xj_to_ch = $(m.Xj_to_ch)")
    println(io, "qJ = $(m.qJ)")
    println(io, "vJ = $(m.vJ)")
    println(io, "cJ = $(m.cJ)")
end

#-------------------------------------------------------------------------------
mutable struct System

end


#-------------------------------------------------------------------------------
# add a single body
function AddBody(id::Int, cf::ConfigBody, b::SingleBody)
    # hierarchy info
    b.bid = id

    # verts
    b.nverts = cf.nverts
    b.verts = zeros(Float64, b.nverts, 3)
    b.verts[:,1] = cf.verts[:,2]
    b.verts[:,3] = cf.verts[:,1]
    b.verts_i = b.verts

    # x_c
    Xc = 0.; Zc = 0.; A = 0.; Ix = 0.; Iz = 0.; Ixz = 0.
    for i = 1:b.nverts
        xj = b.verts[i,1]
        zj = b.verts[i,3]
        if i == b.nverts
            xjp1 = b.verts[1,1]
            zjp1 = b.verts[1,3]
        else
            xjp1 = b.verts[i+1,1]
            zjp1 = b.verts[i+1,3]
        end
        fact = zj*xjp1 - zjp1*xj
        Xc = Xc + (xj + xjp1)*fact
        Zc = Zc + (zj + zjp1)*fact
        A  = A  + 0.5*fact
    end
    Xc = Xc/(6.*A); Zc = Zc/(6.*A)
    b.x_c = [Xc, 0., Zc]

    # mass in scalar, inertia_c in 6d form
    b.mass  = cf.ρ*A
    for i = 1:b.nverts
        # make (Zc,Xc) the origin
        xj = b.verts[i,1] - Xc
        zj = b.verts[i,3] - Zc
        if i == b.nverts
            xjp1 = b.verts[1,1] - Xc
            zjp1 = b.verts[1,3] - Zc
        else
            xjp1 = b.verts[i+1,1] - Xc
            zjp1 = b.verts[i+1,3] - Zc
        end
        fact = zj*xjp1 - zjp1*xj
        Ix =   Ix + (zj^2+zj*zjp1+zjp1^2)*fact
        Iz =   Iz + (xj^2+xj*xjp1+xjp1^2)*fact
        Ixz = Ixz + (xj*zjp1+2.*xj*zj+2.*xjp1*zjp1+xjp1*zj)*fact
    end
    Ix = Ix/12.
    Iz = Iz/12.
    Iy = Ix + Iz
    Ixz = Ixz/24.
    inertia_3d = cf.ρ*[Ix 0. -Ixz; 0. Iy 0.; -Ixz 0. Iz]
    mass_3d = b.mass*eye(3)
    b.inertia_c = [inertia_3d zeros(Float64,3,3);
                   zeros(Float64,3,3) mass_3d]

    # Xb_to_c
    b.Xb_to_c = TransMatrix([zeros(Float64,3);b.x_c])

    # inertia_i in inertial frame
    b.inertia_i = b.Xb_to_c'*b.inertia_c*b.Xb_to_c

    # init q, v, c
    b.q = zeros(Float64,6)
    b.v = zeros(Float64,6)
    b.c = zeros(Float64,6)

    return b
end

#-------------------------------------------------------------------------------
# add a single joint
function AddJoint(cf::ConfigJoint, j::SingleJoint)
    return j
end

function AssembleSystem()

end






end
