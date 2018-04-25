module Utils

export MassCenter

using ..ConstructSystem

function MassCenter(bs::Vector{SingleBody}, sys::System)
    center = zeros(Float64, 3)
    for i = 1:sys.nbody
        center += 0.5*(bs[i].verts_i[2,:] + bs[i].verts_i[3,:])
    end
    center /= sys.nbody
    return center
end












end
