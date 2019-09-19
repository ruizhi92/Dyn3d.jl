# RigidBodyDynamics contain HERK matrix blocks and right hand side
module RigidBodyDynamics

export HERKFuncM, HERKFuncGT, HERKFuncG, HERKFuncf, HERKFuncgti

using Dyn3d
using LinearAlgebra

"""
    HERKFuncM(sys)

HERKFuncM constructs the input function M for HERK method. It returns the
collected inertia matrix of all body in each of their own body coord
"""
function HERKFuncM(sys::System)
    return sys.Ib_total
end


"""
    HERKFuncf(bs, js, sys, f_exi; ρb)

HERKFuncf construct the input function f for HERK method.
It returns a forcing term, which is a summation of bias force term and
joint spring-damper forcing term. The bias term includes the change of
inertia effect, together with gravity and external force. If the body is
a 2d body (for example a cylinder) in the fluid, buoyancy is also taken into account.
"""
function HERKFuncf(bs::Vector{SingleBody}, js::Vector{SingleJoint}, sys::System,
    f_exi::Array{Float64,2}; influid::Bool=false, bodydim::Int=1, ρb::Float64=0.0)

    # pointer to pre-allocated array
    @getfield sys.pre_array (p_total, τ_total, p_bias, f_g, f_ex, r_temp,
        Xic_to_i, A_total, la_tmp1, la_tmp2)

    # determine buoyancy, with consideration of
    if influid == false || bodydim == 1
        ρ_reduced = 1.0
    else
        ρb == 0.0 ? ρ_reduced = 0.0 : ρ_reduced = 1-1/ρb
    end

    # compute bias force, gravity and external force
    for i = 1:sys.nbody
        # bias force
        p_bias = Mfcross(bs[i].v, (bs[i].inertia_b*bs[i].v))
        # gravity in inertial center coord
        f_g = ρ_reduced*bs[i].mass*[zeros(Float64, 3); sys.g]
        # get transform matrix from x_c in inertial frame to the origin of
        # inertial frame
        r_temp = [zeros(Float64, 3); -bs[i].x_c]
        r_temp = bs[i].Xb_to_i*r_temp
        r_temp = [zeros(Float64, 3); -bs[i].x_i + r_temp[4:6]]
        Xic_to_i = TransMatrix(r_temp,la_tmp1,la_tmp2)
        # transform gravity force
        f_g = bs[i].Xb_to_i'*inv(Xic_to_i')*f_g
        # input external force f_exi described in inertial coord
        f_ex = bs[i].Xb_to_i'*f_exi[i,:]

        # add up
        p_total[6i-5:6i] = p_bias - (f_g + f_ex)
    end

    # construct τ_total, this is related only to spring force.
    # τ is only determined by whether the dof has resistance(damp and
    # stiff) or not. Both active dof and passive dof can have τ term
    for i = 1:sys.nbody, k = 1:js[i].nudof
        # find index of the dof in the unconstrained list of this joint
        dofid = js[i].joint_dof[k].dof_id
        τ_total[js[i].udofmap[k]] = -js[i].joint_dof[k].stiff*js[i].qJ[dofid] -
                                    js[i].joint_dof[k].damp*js[i].vJ[dofid]
    end

    # construct A_total to take in parent-child hierarchy
    for i = 1:sys.nbody
        # fill in parent joint blocks
        A_total[6i-5:6i, 6i-5:6i] = Matrix{Float64}(I, 6, 6)
        # fill in child joint blocks except for those body whose nchild=0
        for child_count = 1:bs[i].nchild
            chid = bs[i].chid[child_count]
            A_total[6i-5:6i, 6chid-5:6chid] = - bs[chid].Xp_to_b
        end
    end

    # collect all together
    return -p_total + A_total*sys.S_total*τ_total
end


"""
    HERKFuncGT(bs, sys)

HERKFuncGT constructs the input function GT for HERK method. It returns the
force constraint matrix acting on Lagrange multipliers, including the information
of both parent-child hierarchy and frame transform between adjacent bodies.
"""
function HERKFuncGT(bs::Vector{SingleBody}, sys::System)

    # pointer to pre-allocated array
    @getfield sys.pre_array (A_total,)

    # construct A_total to take in parent-child hierarchy
    for i = 1:sys.nbody
        # fill in parent joint blocks
        A_total[6i-5:6i, 6i-5:6i] = Matrix{Float64}(I, 6, 6)
        # fill in child joint blocks except for those body whose nchild=0
        for child_count = 1:bs[i].nchild
            chid = bs[i].chid[child_count]
            A_total[6i-5:6i, 6chid-5:6chid] = - (bs[chid].Xp_to_b)'
        end
    end
    return A_total*sys.T_total
end


"""
    HERKFuncG(bs, sys)

HERKFuncG constructs the input function G for HERK method.
It returns the motion constraint matrix acting on all body's velocity.
These constraints arise from body velocity relation in each body's local
body coord, for example if body 2 and 3 are connected then:
   v(3) = vJ(3) + X2_to_3*v(2)
"""
function HERKFuncG(bs::Vector{SingleBody}, sys::System)

    # pointer to pre-allocated array
    @getfield sys.pre_array (B_total,)

    # construct B_total to take in parent-child hierarchy
    for i = 1:sys.nbody
        # fill in child body blocks
        B_total[6i-5:6i, 6i-5:6i] = Matrix{Float64}(I, 6, 6)
        # fill in parent body blocks except for those body whose pid=0
        if bs[i].pid != 0
            pid = bs[i].pid
            B_total[6i-5:6i, 6pid-5:6pid] = - bs[i].Xp_to_b
        end
    end
    return (sys.T_total')*B_total
end


"""
    HERKFuncgti(js, sys)

HERKFuncgti returns all the collected prescribed active velocity of joints
at given time.
"""
function HERKFuncgti(js::Vector{SingleJoint}, sys::System, t::T) where T <: AbstractFloat

    # pointer to pre-allocated array
    @getfield sys.pre_array (v_gti, va_gti)

    # give actual numbers from calling motion(t)
    for i = 1:sys.na
        jid = sys.kinmap[i,1]
        dofid = sys.kinmap[i,2]
        _, va_gti[i] = js[jid].joint_dof[dofid].motion(t)
    end

    v_gti[sys.udof_a] = va_gti
    return -(sys.T_total')*v_gti
end


end
