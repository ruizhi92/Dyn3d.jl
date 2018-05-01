module UpdateSystem

export  UpdatePosition!, UpdateVelocity!, InitSystem!

# use registered packages
using DocStringExtensions

# abstract type Velocity end
# abstract type Position end

# import self-defined modules
using ..ConstructSystem
using ..SpatialAlgebra

#-------------------------------------------------------------------------------
function Jcalc(kind::String, qJ::Vector{T}) where T <: AbstractFloat
    if kind == "revolute" || kind == "prismatic" || kind == "cylindrical" ||
       kind == "planar" || kind == "extended_revolute"
       Xj = TransMatrix(qJ)
   elseif kind == "helical"
       # set fixed screw parameter h
       h = 0.1
       r = h*qJ[4:6]
       θ = qJ[1:3]
       Xj = TransMatrix([θ; r])
   elseif kind == "spherical"
       # use quartonians
       error("Under construction")
   elseif kind == "free"
       error("Under construction")
   end
   return Xj
end

#-------------------------------------------------------------------------------
function UpdatePosition!(bs::Vector{SingleBody}, js::Vector{SingleJoint},
    sys::System, qJ::Vector{T} = [2018.]) where T <: AbstractFloat
"""
    1. In InitSystem, apply the part which is the original embed_system
    2. In the middle of HERK, updated qJ in js and update related quantities in
    bs and js.
"""
    if qJ != convert(Vector{T},[2018])
        count = 0
        for i = 1:sys.njoint
            js[i].qJ = qJ[count+1:count+6]
            count += 6
        end
    end
"""
    The following is the original EmbedSystem function
    1. Update body chain in the inertial system, including
       Xb_to_i, verts_i and x_i in body system.
    2. Update js[i].Xj by calling Jcalc
    3. Update bs[i].Xp_to_b using Xp_to_b = Xj_to_ch*Xj*Xp_to_j
"""
    #-------------------------------------------------
    # First body
    #-------------------------------------------------
    # starting from joint 1, which is connected to the inertial system
    # (always true) and is the parent of every other body.
    # js[1].Xj
    js[1].Xj = Jcalc(js[1].joint_type, js[1].qJ)
    # bs[1].Xp_to_b and bs[1].Xb_to_i
    bs[1].Xp_to_b = js[1].Xj_to_ch*js[1].Xj*js[1].Xp_to_j
    bs[1].Xb_to_i = inv(bs[1].Xp_to_b)
    # bs[1].x_i
    q_temp = zeros(Float64, 6)
    q_temp[4:6] = js[1].qJ[4:6]
    q_temp = inv(js[1].Xp_to_j)*q_temp
    bs[1].x_i = js[1].shape1[4:6] + q_temp[4:6]
    #-------------------------------------------------
    # First to last body following parent-child hierarchy
    #-------------------------------------------------
    for i = 1:sys.njoint
        # update verts_i
        for k = 1:bs[i].nverts
            q_temp = zeros(Float64, 6)
            q_temp[4:6] = bs[i].verts[k,:]
            q_temp = bs[i].Xb_to_i*q_temp
            bs[i].verts_i[k,:] = q_temp[4:6] + bs[i].x_i
        end
        # for this joint, loop through every child of it.
        for k = 1:bs[i].nchild
            chid = bs[i].chid[k]
            # update js[i].Xj
            js[chid].Xj = Jcalc(js[chid].joint_type, js[chid].qJ)
            # bs[chid].Xp_to_b
            bs[chid].Xp_to_b = js[chid].Xj_to_ch*js[chid].Xj*
                               js[chid].Xp_to_j
            # bs[chid].Xb_to_i
            bs[chid].Xb_to_i = bs[i].Xb_to_i*inv(bs[chid].Xp_to_b)
            # update x_0 for this child in the inertial system
            # step 1: find the vector to account for shape1(shape1 is expressed
            #         in the parent joint coord)
            q_temp = [ zeros(Float64, 3); js[chid].shape1[4:6] ]
            q_temp = bs[i].Xb_to_i*q_temp
            x_temp = q_temp[4:6] + bs[i].x_i
            # step 2: find the vector to account for joint displacement
            #         (qj is expressed in the child joint coord)
            q_temp = [ zeros(Float64, 3); js[chid].qJ[4:6] ]
            q_temp = bs[i].Xb_to_i*inv(js[chid].Xp_to_j)*q_temp
            x_temp = x_temp + q_temp[4:6]
            # step 3: find the vector to accout for shape2(shape2 is expressed
            #         in the child joint coord)
            q_temp = [ zeros(Float64, 3); -js[chid].shape2[4:6] ]
            q_temp = bs[chid].Xb_to_i*q_temp
            x_temp = x_temp + q_temp[4:6]
            # assign to bs[chid].x_i
            bs[chid].x_i = x_temp
        end
    end
    return bs, js, sys
end

#-------------------------------------------------------------------------------
function UpdateVelocity!(bs::Vector{SingleBody}, js::Vector{SingleJoint},
    sys::System, v::Vector{T}) where T <: AbstractFloat
"""
    With updated body velocity in the middle of HERK, update bs and js,
    return joint velocity in one array.
"""
    # update bs[i].v using input argument v
    count = 0
    for i = 1:sys.nbody
        bs[i].v = v[count+1:count+6]; count += 6
    end
    # update js[i].vJ
    for i = 1:sys.njoint
        pid = bs[i].pid
        if pid != 0
            js[i].vJ = bs[i].v - bs[i].Xp_to_b*bs[pid].v
        else
            js[i].vJ = bs[i].v
        end
        # for planar type joints, we rotate the joint velocity from Fs back to
        # Fp so that the integrated result in q described in Fp coord,
        # which can be used directly as transformation matrix.
        if js[i].joint_type == "planar" ||
            js[i].joint_type == "extended_revolute"
            rot = zeros(T, 6, 6)
            rot[1:3, 1:3] = js[i].Xj[1:3, 1:3]
            rot[4:6, 4:6] = rot[1:3, 1:3]
            js[i].vJ = (rot')*js[i].vJ
        end
    end
    # return joint velocity in an array
    count = 0; vJ = zeros(T, size(v))
    for i = 1:sys.njoint
        vJ[count+1:count+6] = js[i].vJ; count += 6
    end
    return bs, js, sys, vJ
end

#-------------------------------------------------------------------------------
function InitSystem!(bs::Vector{SingleBody}, js::Vector{SingleJoint},
    sys::System)
"""
    InitSystem initialize the joint-body chain by assining values to
    bs[i].v and js[i].qJ at time=0. Body velocity are calculated from joint
    velocities, which are got through articulated body method, zero-out
    total initial momentum.
"""
    # insert active motion
    for i = 1:sys.njoint, k = 1:js[i].na
        act= js[i].udof_a[k]
        idx = js[i].i_udof_a[k]
        js[i].qJ[act] = js[i].joint_dof[idx].motion(0.0)[1]
    end
    # update body chain position
    bs, js, sys = UpdatePosition!(bs, js, sys)
    # zero-out total initial momentum by assigning extra velocity to joint 1's
    # passive dof.
    # pass 1, from body 1 to body n
    for i = 1:sys.nbody
        # initialize the articulated inertia of each body to be equal to its own
        bs[i].Ib_A = bs[i].inertia_b
        # initialize the joint momentum term
        bs[i].pA = zeros(Float64, 6)
    end
    # pass 2, from body n to body 1
    # If the parent is not the base, add the composite inertia of this body
    # (in the coordinate system of the parent) to the inertia of its parent
    for i = sys.nbody:-1:1
        pid = bs[i].pid
        if pid != 0
            Xp_to_b = bs[i].Xp_to_b
            Ib_A_rest = bs[i].Ib_A
            pA_rest = bs[i].pA + Ib_A_rest*js[i].vJ
            bs[pid].Ib_A += (Xp_to_b')*Ib_A_rest*Xp_to_b
            bs[pid].pA += (Xp_to_b')*pA_rest
        end
    end
    # pass 3, from body 1 to body n
    # compute the velocity of passive degrees of freedom in joint from inertial
    # system to body 1, by zeroing the overall system momentum.
    if js[1].np > 0
        Pup = js[1].S[:, js[1].i_udof_p]
        Ptemp = (Pup')*bs[1].Ib_A*Pup
        udof_p = js[1].udof_p
        js[1].vJ[udof_p] = -inv(Ptemp)*(Pup')*bs[1].pA
    end
    # get body[i].v from updated js[i].vJ
    for i = 1:sys.nbody
        pid = bs[i].pid
        if pid != 0
            bs[i].v = js[i].vJ + bs[i].Xp_to_b*bs[pid].v
        else
            bs[i].v = js[i].vJ
        end
    end
    # construct first solution
    qJ_total = zeros(Float64,6*sys.njoint)
    v_total = zeros(Float64,6*sys.nbody)
    for i = 1:sys.nbody
        qJ_total[6i-5:6i] = js[i].qJ
        v_total[6i-5:6i] = bs[i].v
    end
    soln = Soln(0.0, sys.num_params.dt, qJ_total, v_total)

    return bs, js, sys, soln
end







end
