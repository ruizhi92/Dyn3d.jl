module JointType

export ChooseJoint#, Prismatic, Cylindrical, Planar, Spherical, Free

mutable struct ChooseJoint
    nudof::Int
    ncdof::Int
    udof::Union{Vector{Int},Int}
    cdof::Union{Vector{Int},Int,Nullable{Int}}
    S::Union{Vector{Int},Array{Int,2}}
    T::Union{Vector{Int},Array{Int,2},Nullable{Int}}
end

function ChooseJoint(kind="revolute")
    nudof = 1
    ncdof = 5
    udof = 3
    cdof = [1, 2, 4, 5, 6]
    S = eye(Int,6)[:,udof]
    T = eye(Int,6)[:,cdof]
    return Joint(nudof, ncdof, udof, cdof, S, T)
end

function ChooseJoint(kind="prismatic")
    nudof = 1
    ncdof = 5
    udof = 6
    cdof = [1, 2, 3, 4, 5]
    S = eye(Int,6)[:,udof]
    T = eye(Int,6)[:,cdof]
    return Joint(nudof, ncdof, udof, cdof, S, T)
end

function ChooseJoint(kind="cylindrical")
    nudof = 2
    ncdof = 4
    udof = [3, 6]
    cdof = [1, 2, 4, 5]
    S = eye(Int,6)[:,udof]
    T = eye(Int,6)[:,cdof]
    return Joint(nudof, ncdof, udof, cdof, S, T)
end

function ChooseJoint(kind="planar")
    nudof = 3
    ncdof = 3
    udof = [3, 4, 5]
    cdof = [1, 2, 6]
    S = eye(Int,6)[:,udof]
    T = eye(Int,6)[:,cdof]
    return Joint(nudof, ncdof, udof, cdof, S, T)
end

function ChooseJoint(kind="spherical")
    nudof = 3
    ncdof = 3
    udof = [1, 2, 3]
    cdof = [4, 5, 6]
    S = eye(Int,6)[:,udof]
    T = eye(Int,6)[:,cdof]
    return Joint(nudof, ncdof, udof, cdof, S, T)
end

function ChooseJoint(kind="free")
    nudof = 6
    ncdof = 0
    udof = [1, 2, 3, 4, 5, 6]
    cdof = Nullable{Int}()
    S = eye(Int,6)[:,udof]
    T = Nullable{Int}()
    return Joint(nudof, ncdof, udof, cdof, S, T)
end

a = ChooseJoint("prismatic")
a = ChooseJoint("free")
