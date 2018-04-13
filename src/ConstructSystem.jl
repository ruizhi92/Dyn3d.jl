module ConstructSystem

# export
export SingleBody, SingleJoint, AddBody, AddJoint, AssembleSystem

# use registered packages
using DocStringExtensions
import Base: show

# import self-defined modules
include("JointType.jl")
using .JointType
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
    udof::Vector{Int} # newline in constructor
    cdof::Union{Vector{Int},Nullable{Int}}
    udof_p::Union{Vector{Int},Nullable{Int}} # newline in constructor
    udof_a::Union{Vector{Int},Nullable{Int}}
    i_udof_p::Union{Vector{Int},Nullable{Int}} # newline in constructor
    i_udof_a::Union{Vector{Int},Nullable{Int}}
    udof_map::Vector{Int}
    # dof info modified by HERK
    nudof_HERK::Int
    ncdof_HERK::Int
    udof_HERK::Union{Vector{Int},Nullable{Int}}
    cdof_HERK::Union{Vector{Int},Nullable{Int}}
    # joint basis matrix
    S::Union{Vector{Int},Array{Int,2}}
    T::Union{Vector{Int},Array{Int,2},Nullable{Int}}
    T_HERK::Union{Vector{Int},Array{Int,2},Nullable{Int}}
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
    Vector{Int}(0),Vector{Int}(0),
    Vector{Int}(0),Vector{Int}(0),
    Vector{Int}(0),Vector{Int}(0),Vector{Int}(0),
    0,0,Vector{Int}(0),Vector{Int}(0),
    Array{Int,2}(0,0),Array{Int,2}(0,0),Array{Int,2}(0,0),
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
function AddJoint(id::Int, cf::ConfigJoint, j::SingleJoint)
    # hierarchy info
    j.jid = id
    j.joint_type = cf.joint_type
    j.body1 = cf.body1
    j.shape1 = cf.shape1
    j.shape2 = cf.shape2
    # dof info
    j.joint_dof = cf.joint_dof
    # fill in nudof, ncdof, udof, cdof, S, T
    choosen = ChooseJoint(j.joint_type)
    j.nudof = choosen.nudof
    j.ncdof = choosen.ncdof
    j.udof = choosen.udof
    j.cdof = choosen.cdof
    j.S = choosen.S
    j.T = choosen.T
    # np and na
    j.np = 0; j.na = 0
    for i = 1:j.nudof
        if j.joint_dof[i].dof_type == "passive" j.np += 1
        elseif j.joint_dof[i].dof_type == "active" j.na += 1
        else error("dof_type has to be active or passive") end
    end
    # udof_p and i_udof_p
    count = 1
    if j.np != 0
        j.udof_p = Vector{Int}(j.np)
        j.i_udof_p = Vector{Int}(j.np)
        for i = 1:j.nudof
            if j.joint_dof[i].dof_type == "passive"
                j.udof_p[count] = j.joint_dof[i].dof_id
                j.i_udof_p[count] = i
                count += 1
            end
        end
    end
    # udof_a and i_udof_a
    count = 1
    if j.na != 0
        j.udof_a = Vector{Int}(j.na)
        j.i_udof_a = Vector{Int}(j.na)
        for i = 1:j.nudof
            if j.joint_dof[i].dof_type == "active"
                j.udof_a[count] = j.joint_dof[i].dof_id
                j.i_udof_a[count] = i
                count += 1
            end
        end
    end
    # nudof_HERK and ncdof_HERK
    j.nudof_HERK = j.nudof; j.ncdof_HERK = j.ncdof
    for i = 1:j.nudof
        if j.joint_dof[i].dof_type == "active"
            j.nudof_HERK -= 1; j.ncdof_HERK += 1
        end
    end
    # udof_HERK, modified by active motion
    if j.nudof_HERK != 0
        j.udof_HERK = Vector{Int}(j.nudof_HERK)
        count = 1
        for i = 1:6
            if j.ncdof != 0
                if !(i in j.cdof) j.udof_HERK[count] = i; count += 1 end
            elseif j.na != 0
                if !(i in j.udof_a) j.udof_HERK[count] = i; count += 1 end
            end
        end
    end
    # cdof_HERK, modified by active motion
    if j.ncdof_HERK != 0
        j.cdof_HERK = Vector{Int}(j.ncdof_HERK)
        count = 1
        for i = 1:6
            if j.ncdof != 0
                if (i in j.cdof) j.cdof_HERK[count] = i; count += 1 end
            end
            if j.na != 0
                if (i in j.udof_a) j.cdof_HERK[count] = i; count += 1 end
            end
        end
    end
    # T_HERK
    if j.ncdof_HERK != 0
        j.T_HERK = zeros(Int,6,j.ncdof_HERK)
        count = 1
        for i = 1:j.ncdof_HERK j.T_HERK[j.cdof_HERK[i],i] = 1 end
    end
    # Xp_to_j, Xj_to_ch
    j.Xp_to_j = TransMatrix(j.shape1)
    j.Xj_to_ch = TransMatrix(j.shape2)
    # qJ, vJ and cJ
    j.qJ = cf.qJ_init
    j.vJ = zeros(Float64,6)
    j.cJ = zeros(Float64,6)

    return j
end

function AssembleSystem()

end






end
