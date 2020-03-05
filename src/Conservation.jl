# Currently only support 1d body moving in 2d space[x,y]
module Conservation

using LinearAlgebra
using Dyn3d

export momentum!, force!, energy!, work!

"""
    momentum!(dict,bd)

- `dict` : structure of terms in the energy conservation. Includes keys of
     "mlag", "mlag_temp", "mgra", "mgra_temp", "mkin"
- `bd` : current BodyDyn structure

Computes total body's momentum in inertia frame. This function only needs to be
called once at the end stage in RK scheme.
"""
function momentum!(dict::Dict,bd::BodyDyn)

    @getfield bd (bs,js,sys)

    # buffer
    kin = zeros(2)

    # --------- body momentum ---------
    for i = 1:sys.nbody
        kin .+= (inv(bs[i].Xb_to_i)'*(bs[i].inertia_b*bs[i].v))[4:5]
    end

    # output
    dict["mkin"][end] .= kin
    return dict
end

"""
    force!(dict,bd,λ,k)

- `dict` : structure of terms in the energy conservation. Includes keys of
    "mlag", "mlag_temp", "mgra", "mgra_temp", "mkin"
- `bd` : current BodyDyn structure
- `λ` : force by body on the constrained dof of joints
- `k` : k-1 represents the current the stage in rk scheme

Computes the force by lagrangian multipliers in inertia frame at one stage in a
timestep of rk scheme if the body-joint system is under active prescribed motion. Also
computes gravity force at one stage. This function needs to be called at every
stage in RK scheme.
"""
function force!(dict::Dict,bd::BodyDyn,λ::Vector{Float64},k::Int)

    @getfield bd (bs,js,sys)

    # parameters
    g = sys.g
    nbody = sys.nbody

    # buffer
    λ6d = zeros(6)
    lag_temp = zeros(2)
    gra_temp = zeros(2)

    # ---------  Lagrangian multipliers ---------
    # flag case, where the bodies are pinned at first hinge and passive for the rest
    if js[1].joint_type == "revolute"
        λ6d[1:2] .= λ[1:2]
        λ6d[4:6] .= λ[3:5]
        lag_temp .+= (inv(bs[1].Xb_to_i)'*λ6d)[4:5]
    end

    # --------- gravity ---------
    if g != zeros(Float64,3)
        for i = 1:nbody
            gra_temp[2] -= bs[i].mass*sys.g[2]
        end
    end

    # output
    dict["mlag_temp"][end][k-1,:] .= lag_temp
    dict["mgra_temp"][end][k-1,:] .= gra_temp
    return dict
end


"""
    energy!(dict,bd,k)

- `dict` : structure of terms in the energy conservation. Includes keys of
    "edam", "edam_temp", "elag", "elag_temp", "ekin", "espr", "egra"
- `bd` : current BodyDyn structure

For the current joint-body connection, computes kinetic energy, spring potential energy
and gravitational potential energy. This function only needs to be called once
at the end stage in RK scheme.
"""
function energy!(dict::Dict,bd::BodyDyn)

    @getfield bd (bs,js,sys)

    # parameters
    g = sys.g
    nbody = sys.nbody

    # buffer
    kin = 0.0
    spr = 0.0
    gra = 0.0

    # --------- body kinetic energy ---------
    for i = 1:nbody
        vtmp = bs[i].v
        kin += 0.5*(bs[i].inertia_b*vtmp)'*vtmp
    end

    # --------- spring potential energy ---------
    # This is valid for:
    # 1. passive dof
    # 2. active dof with non-zero stiffness
    for i = 1:nbody
        for m = 1:js[i].nudof
            dofid = js[i].joint_dof[m].dof_id
            spr += 0.5*js[i].joint_dof[m].stiff*(js[i].qJ[dofid]^2)
        end
    end

    # --------- gravitational potential energy ---------
    if g != zeros(Float64,3)
        for i = 1:nbody
            # y0 = 0.5*(bdhist[1].bs[i].verts_i[2,2] + bdhist[1].bs[i].verts_i[3,2])
            y = 0.5*(bs[i].verts_i[2,2] + bs[i].verts_i[3,2])
            gra += bs[i].mass*abs(g[2])*(y)
        end
    end

    # output
    dict["ekin"][end] = kin
    dict["espr"][end] = spr
    dict["egra"][end] = gra
    return dict
end

"""
    work!(dict,bd,λ,k)

- `dict` : structure of terms in the energy conservation. Includes keys of
    "edam", "edam_temp", "elag", "elag_temp", "ekin", "espr", "egra"
- `bd` : current BodyDyn structure
- `λ` : force by body on the constrained dof of joints
- `k` : k-1 represents the current the stage in rk scheme

Computes the work done by lagrangian multipliers within one stage in a timestep
of rk scheme if the body-joint system is under active prescribed motion. Also
computes the energy extracted by damper within one stage. This function needs to
be called at every stage in RK scheme.
"""
function work!(dict::Dict,bd::BodyDyn,λ::Vector{Float64},k::Int)

    @getfield bd (bs,js,sys)

    # parameters
    nbody = sys.nbody

    # buffer
    dam_temp = 0.0
    lag_temp = 0.0

    # --------- power extracted by damper ---------
    for i = 1:nbody
        vJtmp = js[i].vJ
        for m = 1:js[i].nudof
            dofid = js[i].joint_dof[m].dof_id
            dam_temp += js[i].joint_dof[m].damp*(vJtmp[dofid]^2)
        end
    end

    # --------- work by Lagrangian force to enforce active motion ---------
    atmp = 0 # temporary active dof
    for i = 1:nbody
        na = js[i].na
        if na > 0
            udof_a = js[i].udof_a
            anker = atmp
            for j = 1:na
                anker = atmp + udof_a[j]
                lag_temp += λ[anker]*js[i].vJ[udof_a[j]]
            end
        end
        atmp += 6-js[i].np
    end

    # output
    dict["edam_temp"][end][k-1] = dam_temp
    dict["elag_temp"][end][k-1] = lag_temp
    return dict
end

end
