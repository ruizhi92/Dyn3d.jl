module Utils

export MassCenter, VertsHistory

using ..ConstructSystem

#-------------------------------------------------------------------------------
function MassCenter(bs::Vector{SingleBody}, sys::System)
    center = zeros(Float64, 3)
    for i = 1:sys.nbody
        center += 0.5*(bs[i].verts_i[2,:] + bs[i].verts_i[3,:])
    end
    center /= sys.nbody
    return center
end

#-------------------------------------------------------------------------------
function VertsHistory(nbody::Int, bs::Vector{SingleBody})
    verts_i = Array{Float64}(nbody,bs[1].nverts,3)
    for i = 1:nbody
        verts_i[i,:,:] = bs[i].verts_i
    end
    return verts_i
end



end
