module ConstructSystem

# export
export SingleBody, SingleJoint, System, Soln,
       AddBody, AddJoint, AssembleSystem!

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
    3. AssembleSystem!
"""

#-------------------------------------------------------------------------------
mutable struct Soln{T}
    # current time and the next timestep
    t::T
    dt::T
    # position, velocity, acceleration and Lagrange multipliers
    qJ::Vector{T}
    v::Vector{T}
    v̇::Vector{T}
    λ::Vector{T}
end

Soln(t, dt, q, v) = Soln(t, dt, q, v,
    Vector{typeof(t)}(0), Vector{typeof(t)}(0))

Soln(t) = Soln(t, 0., Vector{typeof(t)}(0), Vector{typeof(t)}(0),
        Vector{typeof(t)}(0), Vector{typeof(t)}(0))

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
    inertia_b::Array{Float64,2}
    # transform matrix
    Xb_to_c::Array{Float64,2}
    Xb_to_i::Array{Float64,2}
    Xp_to_b::Array{Float64,2}
    # motion vectors
    q::Vector{Float64}
    v::Vector{Float64}
    v̇::Vector{Float64}
    # articulated body info
    pA::Vector{Float64}
    Ib_A::Array{Float64,2}
end

# outer constructor
function SingleBody()
    body = SingleBody(
    0,0,Vector{Int}(0),0,
    0,Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    Vector{Float64}(0),Vector{Float64}(0),
    0,Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    Array{Float64,2}(0,0),Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    Vector{Float64}(0),Vector{Float64}(0),Vector{Float64}(0),
    Vector{Float64}(0),Array{Float64,2}(0,0)
    )
    return body
end

function show(io::IO, ::MIME"text/plain", m::SingleBody)
    println(io, "body_id = $(m.bid)", ", parent_id = $(m.pid)",
            ", nchild = $(m.nchild)", ", chid = $(m.chid)")
    println(io, "nverts = $(m.nverts)", ", verts = $(m.verts)")
    println(io, "verts_i = $(m.verts_i)")
    println(io, "x_c = $(m.x_c)", ", x_i = $(m.x_i)")
    println(io, "mass = $(m.mass)")
    println(io, "inertia_c = $(m.inertia_c)")
    println(io, "inertia_b = $(m.inertia_b)")
    println(io, "Xb_to_c = $(m.Xb_to_c)")
    println(io, "Xb_to_i = $(m.Xb_to_i)")
    println(io, "Xp_to_b = $(m.Xp_to_b)")
    println(io, "q = $(m.q)")
    println(io, "v = $(m.v)")
    println(io, "v̇ = $(m.v̇)")
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
    udofmap::Vector{Int}
    # dof info modified by HERK
    nudof_HERK::Int
    ncdof_HERK::Int
    udof_HERK::Union{Vector{Int},Nullable{Int}}
    cdof_HERK::Union{Vector{Int},Nullable{Int}}
    cdofmap_HERK::Union{Vector{Int},Nullable{Int}}
    # joint basis matrix
    S::Union{Vector{Int},Array{Int,2}}
    T::Union{Vector{Int},Array{Int,2},Nullable{Int}}
    T_HERK::Union{Vector{Int},Array{Int,2},Nullable{Int}}
    # joint dof
    joint_dof::Vector{Dof}
    # motion vectors
    qJ::Vector{Float64}
    vJ::Vector{Float64}
    v̇J::Vector{Float64}
    # transform matrix
    Xj::Array{Float64,2}
    Xp_to_j::Array{Float64,2}
    Xj_to_ch::Array{Float64,2}
end

# outer constructor
function SingleJoint()
    joint = SingleJoint(
    0," ",0,Vector{Float64}(0),Vector{Float64}(0),
    0,0,0,0,
    Vector{Int}(0),Vector{Int}(0),
    Vector{Int}(0),Vector{Int}(0),
    Vector{Int}(0),Vector{Int}(0),Vector{Int}(0),
    0,0,Vector{Int}(0),Vector{Int}(0),Vector{Int}(0),
    Array{Int,2}(0,0),Array{Int,2}(0,0),Array{Int,2}(0,0),
    [Dof()], # dof
    Vector{Float64}(0),Vector{Float64}(0),Vector{Float64}(0),
    Array{Float64,2}(0,0),Array{Float64,2}(0,0),Array{Float64,2}(0,0)
    )
    return joint
end

function show(io::IO, m::SingleJoint)
    println(io, "joint_id = $(m.jid)", ", joint_type = $(m.joint_type)",
            ", body1 = $(m.body1)")
    println(io, "shape1 = $(m.shape1)", ", shape2 = $(m.shape2)")
    println(io, "nudof = $(m.nudof)", ", ncdof = $(m.ncdof)",
            ", np = $(m.np)", ", na = $(m.na)")
    println(io, "udof = $(m.udof)", ", cdof = $(m.cdof)")
    println(io, "udof_p = $(m.udof_p)", ", udof_a = $(m.udof_a)")
    println(io, "i_udof_p = $(m.i_udof_p)", ", i_udof_a = $(m.i_udof_a)")
    println(io, "udofmap = $(m.udofmap)")
    println(io, "cdofmap_HERK = $(m.cdofmap_HERK)")
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
    println(io, "v̇J = $(m.v̇J)")
end

#-------------------------------------------------------------------------------
mutable struct PreArray
    # used in UpdatePosition, UpdateVelocity
    q_temp::Vector{Float64}
    x_temp::Vector{Float64}
    rot::Array{Float64,2}
    # used in HERK!
    qJ::Array{Float64,2}
    vJ::Array{Float64,2}
    v::Array{Float64,2}
    v̇::Array{Float64,2}
    λ::Array{Float64,2}
    v_temp::Vector{Float64}
    lhs::Array{Float64,2}
    rhs::Vector{Float64}
    # used in HERKFuncM
    Mᵢ₋₁::Array{Float64,2}
    # used in HERKFuncf
    p_total::Vector{Float64}
    τ_total::Vector{Float64}
    p_bias::Vector{Float64}
    f_g::Vector{Float64}
    f_ex::Vector{Float64}
    r_temp::Vector{Float64}
    Xic_to_i::Array{Float64,2}
    fᵢ₋₁::Vector{Float64}
    # used in HERKFuncf, HERKFuncGT
    A_total::Array{Float64,2}
    GTᵢ₋₁::Array{Float64,2}
    # used in HERKFuncG
    B_total::Array{Float64,2}
    Gᵢ::Array{Float64,2}
    # used in HERKFuncgti
    v_gti::Vector{Float64}
    va_gti::Vector{Float64}
    gtiᵢ::Vector{Float64}
end

PreArray() = PreArray(
    # UpdatePosition, UpdateVelocity
    Float64[],Float64[],Array{Float64,2}(0,0),
    # HERK!
    Array{Float64,2}(0,0),Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    Array{Float64,2}(0,0),Array{Float64,2}(0,0),Float64[],
    Array{Float64,2}(0,0),Float64[],
    # HERKFuncM
    Array{Float64,2}(0,0),
    # HERKFuncf
    Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],
    Array{Float64,2}(0,0),Float64[],
    # HERKFuncf, HERKFuncGT
    Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    # HERKFuncG
    Array{Float64,2}(0,0),Array{Float64,2}(0,0),
    # HERKFuncgti
    Float64[],Float64[],Float64[]
)
#-------------------------------------------------------------------------------
mutable struct System
    # general info
    ndim::Int
    nbody::Int
    njoint::Int
    # time info
    time::Float64
    dt::Float64
    # hierarchy info
    ndof::Int
    nudof::Int
    ncdof::Int
    np::Int
    na::Int
    udof::Vector{Int}
    udof_a::Vector{Int}
    udof_p::Vector{Int}
    # hierarchy info modified by HERK
    nudof_HERK::Int
    ncdof_HERK::Int
    udof_HERK::Vector{Int}
    # system info in matrix
    S_total::Array{Int,2}
    T_total::Array{Int,2}
    Ib_total::Array{Float64,2}
    # gravity
    g::Vector{Float64}
    # kinematic map
    kinmap::Array{Int,2}
    # numerical parameters
    num_params::NumParams
    # pre-allocation array
    pre_array::PreArray
end

# outer constructor
System(ndim, nbody, njoint, g, num_params) = System(
    ndim,nbody,njoint,0.,0.,
    0,0,0,0,0,
    Vector{Int}(0),Vector{Int}(0),Vector{Int}(0),
    0,0,Vector{Int}(0),
    Array{Int,2}(0,0),Array{Int,2}(0,0),Array{Float64,2}(0,0),
    g,Array{Int,2}(0,0),num_params,
    PreArray()
)

function show(io::IO, m::System)
    println(io, "ndim = $(m.ndim)", ", njoint = $(m.njoint)",
            ", nbody = $(m.nbody)")
    println(io, "ndof = $(m.ndof)", ", nudof = $(m.nudof)",
            ", ncdof = $(m.ncdof)", ", np = $(m.np)", ", na = $(m.na)")
    println(io, "udof = $(m.udof)")
    println(io, "udof_p = $(m.udof_p)")
    println(io, "udof_a = $(m.udof_a)")
    println(io, "nudof_HERK = $(m.nudof_HERK)",", ncdof_HERK = $(m.ncdof_HERK)")
    println(io, "udof_HERK = $(m.udof_HERK)")
    println(io, "gravity = $(m.g)")
    println(io, "kinmap = $(m.kinmap)")
end

#-------------------------------------------------------------------------------
# add a single body
function AddBody(id::Int, cf::ConfigBody)
    # init a single body object
    b = SingleBody()

    # hierarchy info
    b.bid = id

    # verts
    b.nverts = cf.nverts
    b.verts = zeros(Float64, b.nverts, 3)
    b.verts[:,1] = cf.verts[:,2]
    b.verts[:,3] = cf.verts[:,1]
    b.verts_i = copy(b.verts)

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

    # inertia_b in body frame at point b
    b.inertia_b = b.Xb_to_c'*b.inertia_c*b.Xb_to_c

    # init q, v, c
    b.q = zeros(Float64,6)
    b.v = zeros(Float64,6)
    b.v̇ = zeros(Float64,6)

    return b
end

#-------------------------------------------------------------------------------
# add a single joint
function AddJoint(id::Int, cf::ConfigJoint)
    # init SingleJoint
    j = SingleJoint()

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
    # udof_HERK, modified by active motion
    if j.nudof_HERK != 0
        j.udof_HERK = Vector{Int}(j.nudof_HERK)
        count = 1
        for i = 1:6
            if !(i in j.cdof_HERK) j.udof_HERK[count] = i; count += 1 end
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
    j.qJ = zeros(Float64,6)
    j.qJ[j.udof] = cf.qJ_init
    j.vJ = zeros(Float64,6)
    j.v̇J = zeros(Float64,6)

    return j
end

#-------------------------------------------------------------------------------
# connects body and joint system
function AssembleSystem!(bs::Vector{SingleBody}, js::Vector{SingleJoint},
                        sys::System)
    #-------------------------------------------------
    # re-order bs and js to make the array in id-order
    #-------------------------------------------------
    bstemp = deepcopy(bs)
    for i = 1:sys.nbody
        bs[bstemp[i].bid] = bstemp[i]
    end
    jstemp = deepcopy(js)
    for i = 1:sys.njoint
        js[jstemp[i].jid] = jstemp[i]
    end
    #-------------------------------------------------
    # fill in sys structure
    #-------------------------------------------------
    # scalar hierarchy info
    sys.ndof = 6*sys.njoint
    sys.nudof = 0; sys.ncdof = 0; sys.np = 0; sys.na = 0
    sys.nudof_HERK = 0; sys.ncdof_HERK = 0
    for i = 1:sys.njoint
        sys.nudof += js[i].nudof
        sys.ncdof += js[i].ncdof
        sys.np += js[i].np
        sys.na += js[i].na
        sys.nudof_HERK += js[i].nudof_HERK
        sys.ncdof_HERK += js[i].ncdof_HERK
    end
    # udof
    count = 1
    sys.udof = Vector{Int}(sys.nudof)
    for i = 1:sys.njoint, k = 1:js[i].nudof
        sys.udof[count] = 6*(i-1) + js[i].udof[k]
        count += 1
    end
    # udof_a
    count = 1
    sys.udof_a = Vector{Int}(sys.na)
    for i = 1:sys.njoint, k = 1:js[i].na
        sys.udof_a[count] = 6*(i-1) + js[i].udof_a[k]
        count += 1
    end
    # udof_p
    count = 1
    sys.udof_p = Vector{Int}(sys.np)
    for i = 1:sys.njoint, k = 1:js[i].np
        sys.udof_p[count] = 6*(i-1) + js[i].udof_p[k]
        count += 1
    end
    # udof_HERK
    count = 1
    sys.udof_HERK = Vector{Int}(sys.nudof_HERK)
    for i = 1:sys.njoint, k = 1:js[i].nudof_HERK
        sys.udof_HERK[count] = 6*(i-1) + js[i].udof_HERK[k]
        count += 1
    end
    #-------------------------------------------------
    # fill in info in bs and js structure
    #-------------------------------------------------
    # loop through every joint, find nchild for every body, then allocate
    # bs[i].chid. Also assign bs.pid
    for i = 1:sys.nbody
        bs[i].nchild = 0
        ch_cnt = 1
        for k = 1:sys.njoint
            if js[k].body1 == bs[i].bid bs[i].nchild += 1 end
        end
        # if at least exists one child, allocate chid
        if bs[i].nchild != 0
            bs[i].chid = Vector{Int}(bs[i].nchild)
            for k = 1:sys.njoint
                if js[k].body1 == bs[i].bid
                    bs[i].chid[ch_cnt] = js[k].jid
                    ch_cnt += 1
                end
            end
        end
        # bs[i].pid
        bs[i].pid = js[bs[i].bid].body1
    end
    # js[i].udofmap
    last = 0
    for i = 1:sys.njoint
        if js[i].nudof != 0
                js[i].udofmap = Vector{Int}(js[i].nudof)
                js[i].udofmap = last + [k for k=1:js[i].nudof]
                last = js[i].udofmap[js[i].nudof]
        end
    end
    # js[i].cdofmap_HERK
    last = 0
    for i = 1:sys.njoint
        if js[i].ncdof_HERK != 0
            js[i].cdofmap_HERK = Vector{Int}(js[i].ncdof_HERK)
            js[i].cdofmap_HERK = last + [k for k=1:js[i].ncdof_HERK]
            last = js[i].cdofmap_HERK[js[i].ncdof_HERK]
        end
    end
    # kinmap
    sys.kinmap = Array{Int}(sys.na, 2)
    count = 1
    for i = 1:sys.njoint, k = 1:js[i].na
        sys.kinmap[count,1] = i
        sys.kinmap[count,2] = js[i].i_udof_a[k]
        count += 1
    end
    #-------------------------------------------------
    # change body coord origin from body origin to joint
    #-------------------------------------------------
    for i = 1:sys.njoint
        # store old values
        r_old = js[i].shape2[1:3]
        Xj_to_ch_old = js[i].Xj_to_ch
        rot_old = Xj_to_ch_old[1:3,1:3]
        # The body frame origin now coincides with the joint location
        js[i].shape2[1:3] = 0.0
        # the new transform from parent joint to body is the identity
        js[i].Xj_to_ch = eye(Float64,6)
        # update all the vertices
        for k = 1:bs[i].nverts
            r_temp = -r_old + bs[i].verts[k,:]
            verts_temp = (rot_old')*r_temp
            bs[i].verts[k,:] = verts_temp
        end
        # update body center of mass
        r_temp = -r_old + bs[i].x_c
        x_c_temp = (rot_old')*r_temp
        bs[i].x_c = x_c_temp
        # update bs[i].Xb_to_c and bs[i].inertia_b
        bs[i].Xb_to_c = bs[i].Xb_to_c*Xj_to_ch_old
        bs[i].inertia_b = (Xj_to_ch_old')*bs[i].inertia_b*Xj_to_ch_old
        # update Xp_to_j(related to shape1) for child body
        if bs[i].nchild != 0
            for k = 1:bs[i].nchild
                ck = bs[i].chid[k]
                r_temp = -r_old + js[ck].shape1[4:6]
                js[ck].Xp_to_j = js[ck].Xp_to_j*Xj_to_ch_old
            end
        end
    end

    #-------------------------------------------------
    # Create some system matrix used in HERK
    #-------------------------------------------------
    # diagonal block of sys.Ib_total is the inertia of each body inertia_b
    sys.Ib_total = zeros(Float64, sys.ndof, sys.ndof)
    for i = 1:sys.nbody
        sys.Ib_total[6i-5:6i, 6i-5:6i] = bs[i].inertia_b
    end
    # diagonal block of S_total is S of each joint in joint coord
    sys.S_total = zeros(Int, sys.ndof, sys.nudof)
    for i = 1:sys.njoint
        if js[i].nudof != 0
            sys.S_total[6i-5:6i, js[i].udofmap] = js[i].S
        end
    end
    # diagonal block of T_total is T_HERK of each joint in joint coord
    sys.T_total = zeros(Int, sys.ndof, sys.ncdof_HERK)
    for i = 1:sys.njoint
        if js[i].ncdof_HERK != 0
            sys.T_total[6i-5:6i, js[i].cdofmap_HERK] = js[i].T_HERK
        end
    end

    return bs, js, sys
end


end
