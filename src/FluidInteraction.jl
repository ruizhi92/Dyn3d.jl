module FluidInteraction

export BodyGrid, CutOut2d
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
    f_ex3d::Vector{Vector{Float64}}  # the external force on all body points in inertial grid frame
    f_ex6d::Vector{Float64}  # f_ex3d integrated and transformed to body's origin
end
BodyGrid(bid,np,points) = BodyGrid(bid,np,points,[zeros(3) for i=1:np],
    [zeros(3) for i=1:np],[zeros(3) for i=1:np],zeros(6))

#-------------------------------------------------------------------------------
"""
This function need to be called only once after GenerateBodyGrid for 2d case of
flat plats.

Since both for 2d and 3d cases, bodies are constructed by quadrilateral/triangles,
not lines. Thus for 2d cases where only the line matters, we cut out the info on
the other sides of the plate. Note that verts are formulated in clockwise
direction, with the left-bottom corner as origin.
"""
function CutOut2d(bd::BodyDyn,bgs::Vector{BodyGrid})
    if bd.sys.ndim == 2 && bd.bs[1].nverts == 4
        for i = 1:length(bgs)
            nverts = bd.bs[bgs[i].bid].nverts
            cutout = round(Int,(bgs[i].np-1)/nverts)
            bgs[i].np = round((bgs[i].np-1)/4)+1
            bgs[i].points = bgs[i].points[end:-1:end-cutout]
            bgs[i].q_i = bgs[i].q_i[end:-1:end-cutout]
            bgs[i].v_i = bgs[i].v_i[end:-1:end-cutout]
            bgs[i].f_ex3d = bgs[i].f_ex3d[end:-1:end-cutout]
        end
    end
    return bgs
end

#-------------------------------------------------------------------------------
function GenerateBodyGrid(bd::BodyDyn; np=101)
"""
Given BodyDyn structure, where each body only consists of several verts(usually
4 for quadrilateral and 3 for triangle), return the verts position in inertial
frame of given number of points np by interpolation, of all bodies in the system.
"""
    # here we assume the body chain consists of only 1 body, or several bodies
    # of the same shape
    @getfield bd (bs,sys)

    bodygrids = Vector{BodyGrid}(sys.nbody)
    # for cases with only 1 body, which has more than 4 grid points(like a circle)
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
    @getfield bd (bs, sys)

    # the j-th q_i in body points of a body = bs[i].x_i + Xb_to_i*points[j]
    # the j-th v_i in body points of a body is calculated by transferring to
    # a coordinate that sits at the beginning point of the first body but with
    # zero angle.

    for i = 1:length(bgs)
        b = bs[bgs[i].bid]
        if b.bid == 1
            X_ref = TransMatrix([zeros(Float64,3);b.x_i])
        end
        for j = 1:bgs[i].np
            q_temp = [zeros(Float64, 3); bgs[i].points[j]]
            q_temp = [zeros(Float64, 3); b.x_i] + b.Xb_to_i*q_temp
            bgs[i].q_i[j] = q_temp[4:6]
            v_temp = bs[i].v + [zeros(Float64, 3); cross(bs[i].v[1:3],bgs[i].points[j])]
            bgs[i].v_i[j] = (X_ref*b.Xb_to_i*v_temp)[4:6]
        end
    end
    return bgs
end

#-------------------------------------------------------------------------------
function IntegrateBodyGridDynamics(bd::BodyDyn,bgs::Vector{BodyGrid})
"""
Given external 3d linear fluid force f_ex of each body point contained in updated
bgs structure, do intergral to return integrated 6d body force([torque,force])
exerting on the beginning of current body.
"""
    @getfield bd (bs,sys)
    for i = 1:length(bgs)
        b = bs[bgs[i].bid]
        bgs[i].f_ex6d = zeros(6)
        for j = 1:bgs[i].np
            # linear force in inertial grid coord
            f_temp = [zeros(Float64, 3); bgs[i].f_ex3d[j]]
            # get transform matrix from grid points in inertial frame to the origin of inertial frame
            r_temp = [zeros(Float64, 3); -(bgs[i].points[1]-bgs[i].points[j])]
            r_temp = b.Xb_to_i*r_temp
            r_temp = [zeros(Float64, 3); -b.x_i + r_temp[4:6]]
            Xic_to_i = TransMatrix(r_temp)
            # transform force to body local frame
            f_temp = b.Xb_to_i'*inv(Xic_to_i')*f_temp
            # println("j= ",j," points[j]= ",bgs[i].points[j]," f_temp= ",f_temp)
            bgs[i].f_ex6d += f_temp
        end
    end
    return bgs
end


end
