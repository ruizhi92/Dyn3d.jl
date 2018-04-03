module ConfigDataType

export config_body, config_joint

# configuration properties of a single body
struct config_body
    nbody::Int
    nverts::Int
    verts::Array{Float64,2}
    ρ::Float64
end

# configuration properties of a single joint
struct config_joint
    joint_type::String
    joint_id::Int
    shape::Tuple{Vector{Float64},Vector{Float64}}
    body1::Int
    joint_dof::Vector{dof}
end

# single dof information for every dof in a joint
abstract type dof with
    dof_id::Int
    dof_type::String
    stiff::Float64
    damp::Float64
    motion::motion
end

# motion parameters to describe a motion
struct motion

end

end
