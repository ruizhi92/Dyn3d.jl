module FluidInteraction

export BodyGrid
export AcquireBodyGridKinematics, IntegrateBodyGridDynamics, GenerateBodyGrid

using Dyn3d
using Interpolations

#-------------------------------------------------------------------------------
"""
Design this structure to contain info used by fluid-structure interaction
"""
mutable struct BodyGrid
    bid::Int    # body id in the joint-body chain
    np::Int     # number of grid points on this body
    points::Vector{Vector{Float64}}  # the (x,y,z) coord of points in local body frame
    q_i::Vector{Vector{Float64}}  # the position of all body points in inertial frame
    v_i::Vector{Vector{Float64}}  # the velocity of all body points in inertial frame
    f_ex::Vector{Vector{Float64}}  # the external force on all body points in inertial frame
end
BodyGrid(bid,np,points) = BodyGrid(bid,np,points,Vector{Vector{Float64}}(np),
    Vector{Vector{Float64}}(np),Vector{Vector{Float64}}(np))

#-------------------------------------------------------------------------------
function GenerateBodyGrid(bd::BodyDyn; np=101)
"""
Given BodyDyn structure, where each body only consists of several verts(usually
4 for quadrilateral and 3 for triangle), return the verts position in inertial
frame of given number of points np by interpolation, of all bodies in the system.
"""
    # here we assume the body chain consists of only 1 body, or several bodies
    # of the same shape
    @get bd (bs,js,sys)

    bodygrids = Vector{BodyGrid}(sys.nbody)
    if bs[1].nverts != 3 && bs[1].nverts != 4
        a = bs[1].verts
        bodygrids[1] = BodyGrid(1,bs[1].nverts,[a[i,:] for i =1:size(a,1)])
        return bodygrids
    end

    if (np-1) % bd.bs[1].nverts != 0 error("Number of points can't be divided by system.nverts") end

    bodygrids = Vector{BodyGrid}(sys.nbody)
    for i = 1:sys.nbody
        bid = bs[i].bid
        verts_id = linspace(1, np, bs[bid].nverts+1)
        verts = vcat(bs[bid].verts,bs[bid].verts[1,:]')
        it_x = interpolate((verts_id,), verts[:,1], Gridded(Linear()))
        it_y = interpolate((verts_id,), verts[:,2], Gridded(Linear()))
        it_z = interpolate((verts_id,), verts[:,3], Gridded(Linear()))
        grid = [[it_x[j],it_y[j],it_z[j]] for j=1:np]
        bodygrids[i] = BodyGrid(bid,np,grid)
    end
    return bodygrids
end

#-------------------------------------------------------------------------------
function AcquireBodyGridKinematics(bd::BodyDyn, bgs::Vector{BodyGrid})
"""
Given updated bd structure, which contains 3d bs[i].x_i in inertial frame and
6d bs[i].v of each body in the body local frame, return 3d linear q_i and v_i of
each body point in the inertial frame.
"""
    @get bd (bs, js, sys)

    # the j-th q_i in body points of a body = bs[i].x_i + Xb_to_i*points[j]
    # the j-th v_i in body points of a body is calculated by first calculating
    # the linear velocity of point j in the coordinate of body i, then do a
    # Xb_to_i transformation
    for i = 1:length(bgs)
        b = bs[bgs[i].bid]
        for j = 1:bgs[i].np
            q_temp = [zeros(Float64, 3); bgs[i].points[j]]
            q_temp = [zeros(Float64, 3); b.x_i] + b.Xb_to_i*q_temp
            bgs[i].q_i[j] = q_temp[4:6]
            v_temp = bs[i].v + [zeros(Float64, 3); cross(bs[i].v[1:3],bgs[i].points[j])]
            bgs[i].v_i[j] = (b.Xb_to_i*v_temp)[4:6]
        end
    end
    return bgs
end

#-------------------------------------------------------------------------------
function IntegrateBodyGridDynamics(bgs::Vector{BodyGrid})
"""
Given external 3d linear fluid force f_ex of each body point, do intergral to
return integrated 6d body force([torque,force]) exerting on the beginning of
current body
"""
end




end
