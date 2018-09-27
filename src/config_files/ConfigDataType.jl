module ConfigDataType

export ConfigBody, ConfigJoint, ConfigSystem, Dof, Motions, NumParams

import Base: show

#-------------------------------------------------------------------------------
# motion parameters to describe a motion
mutable struct Motions
    motion_type::String
    motion_params::Vector{Float64}
end

Motions() = Motions("", [])

# function-like object, only need velocity using HERK
function (m::Motions)(t)
    if m.motion_type == "hold"
        q = m.motion_params[1]
        v = 0

    elseif m.motion_type == "velocity"
        q = m.motion_params[1] + t*m.motion_params[2]
        v = m.motion_params[2]

    elseif m.motion_type == "oscillatory"
        amp = m.motion_params[1]
        freq = m.motion_params[2]
        phase = m.motion_params[3]
        arg = 2π*freq*t + phase
        q = amp*cos(arg)
        v = -2π*freq*amp*sin(arg)

    elseif m.motion_type == "ramp_1"
        # Eldredge ramp from 2009 AIAA paper
        # parameters are a and t[4], this motion is not periodic
        a = m.motion_params[1]
        tᵣ = m.motion_params[2:5]
        f(t) = cosh(a*(t-tᵣ[1]))
        g(t) = cosh(a*(t-tᵣ[2]))
        u(t) = cosh(a*(t-tᵣ[3]))
        n(t) = cosh(a*(t-tᵣ[4]))

        ḟ(t) = a*sinh(a*(t-tᵣ[1]))
        ġ(t) = a*sinh(a*(t-tᵣ[2]))
        u̇(t) = a*sinh(a*(t-tᵣ[3]))
        ṅ(t) = a*sinh(a*(t-tᵣ[4]))

        q = log(f(t)*n(t)/(g(t)*u(t)))
        v = g(t)*u(t)/(f(t)*n(t))*
              (-f(t)*n(t)*(ġ(t)*u(t)+u̇(t)*g(t))/(g(t).^2*u(t).^2)
                  + (ḟ(t)*n(t)+ṅ(t)*f(t))/(g(t)*u(t)))

    elseif m.motion_type == "ramp_2"
        a = m.motion_params[1]
        q = 0.5*(tanh(a*t) + 1)
        v = 0.5*a*sech(a*t).^2
    else
        error("This motion type does not exist")
    end

    return q, v
end
#-------------------------------------------------------------------------------
# single dof information for every dof in a joint
mutable struct Dof
    dof_id::Int
    dof_type::String
    stiff::Float64
    damp::Float64
    motion::Motions
end

Dof() = Dof(3, "passive", 0.03, 0.01, Motions())
#-------------------------------------------------------------------------------
# configuration properties of a single body
mutable struct ConfigBody
    nbody::Int
    nverts::Int
    verts::Array{Float64,2}
    ρ::Float64
end


"""
the final plotting direction is:
         ^(y)
         |_____>(z)
     (-x)
the coordinate for verts in 2d input is in [z,x]. So y direction only allow
zero-width body
"""
ConfigBody(nbody) = ConfigBody(nbody, 4,
    [0. 0.; 1. 0.; 1. 1./nbody; 0. 1./nbody], 0.01)

function show(io::IO, m::ConfigBody)
    println(io, " nbody = $(m.nbody)")
    println(io, " nverts = $(m.nverts)")
    println(io, " verts = $(m.verts)")
    println(io, " ρ = $(m.ρ)")
end

#-------------------------------------------------------------------------------
# configuration properties of a single joint
mutable struct ConfigJoint
    njoint::Int
    joint_type::String
    shape1::Vector{Float64}
    shape2::Vector{Float64}
    body1::Int
    joint_dof::Vector{Dof}
    qJ_init::Vector{Float64}
end

ConfigJoint(njoint,joint_type) = ConfigJoint(njoint, joint_type,
    [0., 0., 0., 1./njoint, 0., 0.], zeros(Float64,6),
    0, [Dof()], [0.])

function show(io::IO, m::ConfigJoint)
    println(io, " joint type = $(m.joint_type)")
    println(io, " joint position in parent body coord = $(m.shape1)")
    println(io, " joint position in child body coord = $(m.shape2)")
    for i = 1:size(m.joint_dof,1)
        if m.joint_dof[i].dof_type == "passive"
            println(io, " joint unconstrained dof = ",
            "$(m.joint_dof[i].dof_id), under $(m.joint_dof[i].dof_type) motion")
        else
            println(io, " joint unconstrained dof = ",
            "$(m.joint_dof[i].dof_id), under $(m.joint_dof[i].dof_type) ",
            "$(m.joint_dof[i].motion.motion_type) motion")
        end
    end
    println(io, " initial unconstrained dof position = $(m.qJ_init)")
end

#-------------------------------------------------------------------------------
# numerical parameters
mutable struct NumParams
    tf::Float64
    dt::Float64
    scheme::String
    st::Int
    tol::Float64
end

#-------------------------------------------------------------------------------
# config_system parameters
mutable struct ConfigSystem
    ndim::Int
    gravity::Vector{Float64}
    num_params::NumParams
end

end
