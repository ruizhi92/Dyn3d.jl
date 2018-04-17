module JointType

export ChooseJoint

mutable struct ChooseJoint
    nudof::Int
    ncdof::Int
    udof::Vector{Int}
    cdof::Union{Vector{Int},Nullable{Int}}
    S::Union{Vector{Int},Array{Int,2}}
    T::Union{Vector{Int},Array{Int,2},Nullable{Int}}
end

function ChooseJoint(kind)

    if kind == "revolute"
        nudof = 1
        ncdof = 5
        udof = [3]
        cdof = [1, 2, 4, 5, 6]
        S = eye(Int,6)[:,udof]
        T = eye(Int,6)[:,cdof]

    elseif kind == "prismatic"
        nudof = 1
        ncdof = 5
        udof = [6]
        cdof = [1, 2, 3, 4, 5]
        S = eye(Int,6)[:,udof]
        T = eye(Int,6)[:,cdof]

    elseif kind == "cylindrical"
        nudof = 2
        ncdof = 4
        udof = [3, 6]
        cdof = [1, 2, 4, 5]
        S = eye(Int,6)[:,udof]
        T = eye(Int,6)[:,cdof]

    elseif kind == "planar"
        nudof = 3
        ncdof = 3
        udof = [3, 4, 5]
        cdof = [1, 2, 6]
        S = eye(Int,6)[:,udof]
        T = eye(Int,6)[:,cdof]

    elseif kind == "spherical"
        nudof = 3
        ncdof = 3
        udof = [1, 2, 3]
        cdof = [4, 5, 6]
        S = eye(Int,6)[:,udof]
        T = eye(Int,6)[:,cdof]

    elseif kind == "free"
        nudof = 6
        ncdof = 0
        udof = [1, 2, 3, 4, 5, 6]
        cdof = Nullable{Int}()
        S = eye(Int,6)[:,udof]
        T = Nullable{Int}()

    elseif kind == "extended_revolute"
        nudof = 2
        ncdof = 4
        udof = [3, 4]
        cdof = [1, 2, 5, 6]
        S = eye(Int,6)[:,udof]
        T = eye(Int,6)[:,cdof]
    end

    return ChooseJoint(nudof, ncdof, udof, cdof, S, T)
end



end
