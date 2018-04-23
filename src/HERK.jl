module HERK

# export

# use registered packages
using DocStringExtensions

# import self-defined modules
using ..ConstructSystem
using ..SpatialAlgebra

#-------------------------------------------------------------------------------
function HERKScheme(name::String)
"""
HERKScheme provides a set of HERK coefficients, expressed in Butcher table form.
    A: Runge-Kutta matrix in Butcher tableau    c | A
    b: weight vector in Butcher tableau         ------
    c: node vector in Butcher tableau             | b
    s: stage number
    p: method order of accuracy
"""
    # Scheme A of 2-stage HERK in Liska's paper
    if name == "Liska"
        A = [0.0 0.0 0.0;
             0.5 0.0 0.0;
             √3/3 (3.0-√3)/3 0.0]
        b = [(3.0+√3)/6 -√3/3 (3.0+√3)/6]
        c = [0.0, 0.5, 1.0]
        s = 3
        p = 2
    # Brasey-Hairer 3-Stage HERK, table 2
    elseif name == "BH3"
        A = [0.0 0.0 0.0;
             1.0/3 0.0 0.0;
             -1.0 2.0 0.0]
        b = [0.0, 0.75, 0.25]
        c = [0.0, 1.0/3, 1.0]
        s = 3
        p = 3
    # Brasey-Hairer 5-Stage HERK, table 5
    elseif name == "BH5"
        A = [0.0 0.0 0.0 0.0 0.0;
             0.3 0.0 0.0 0.0 0.0;
             (1.0+√6)/30 (11.0-4*√6)/30 0.0 0.0 0.0;
             (-79.0-31*√6)/150 (-1.0-4*√6)/30 (24.0+11*√6)/25 0.0 0.0;
             (14.0+5*√6)/6 (-8.0+7*√6)/6 (-9.0-7*√6)/4 (9.0-√6)/4 0.0]
        b = [0.0, 0.0, (16.0-√6)/36, (16.0+√6)/36, 1.0/9]
        c = [0.0, 0.3, (4.0-√6)/10, (4.0+√6)/10, 1.0]
        s = 5
        p = 4
    else
        error("This HERK scheme doesn't exist now.")
    end
    return A, b, c, s
end

#-------------------------------------------------------------------------------
function HERKMain(tᵢ, qᵢ, vᵢ, δtᵢ, tol, scheme, qᵢ₊₁, vᵢ₊₁, cᵢ₊₁, λᵢ₊₁, δtᵢ₊₁)
"""
    HERKMain is a half-explicit Runge-Kutta solver based on the
    constrained body model in paper of V.Brasey and E.Hairer.
    The following ODE system is being solved:
       | dq/dt = v                     |
       | M(q)*dv/dt = f(q,v) - GT(q)*λ |
       | 0 = G(q)*v + gti(q)           |
    , where GT is similar to the transpose of G.

    Note that this is a index-2 system with two variables q and u.
    Here we denote for a body chain, v stands for body velocity
    and vJ stands for joint velocity. Similarly for q and qJ. λ is the
    constraint on q to be satisfied.

    Specificlly for the dynamics problem, we choose to solve for qJ and v.
    The system of equation is actually:
       | dqJ/dt = vJ                              |
       | M(qJ)*dv/dt = f(qJ,v,vJ) - GT(qJ)*lambda |
       | 0 = G(qJ)*v + gti(qJ)                    |
    So we need a step to calculate v from vJ solved. The motion constraint
    (prescribed active motion) is according to joint, not body.
"""
    return 1



end
#-------------------------------------------------------------------------------
function HERKFuncM(sys::System)
"""
    HERKFuncM constructs the input function M for HERK method.
    It returns the collected inertia matrix of all body in their own body coord
"""
    return sys.Ib_total
end

#-------------------------------------------------------------------------------
function HERKFuncf(bs::Vector{SingleBody}, js::Vector{SingleJoint}, sys::System)
"""
    HERKFuncf construct the input function f for HERK method.
    It returns a forcing term, whchild_count is a summation of bias force term and
    joint spring-damper forcing term. The bias term includes the change of
    inertia effect, together with gravity and external force.
"""
    # allocation
    p_total = zeros(Float64, sys.ndof)
    τ_total = zeros(Float64, sys.nudof)
    A_total = zeros(Float64, sys.ndof, sys.ndof)

    # compute bias force, gravity and external force
    for i = 1:sys.nbody
        # bias force
        p_temp = Mfcross(bs[i].v, (bs[i].inertia_b*bs[i].v))
        # gravity in inertial center coord
        f_g = bs[i].mass*[zeros(Float64, 3); sys.g]
        # get transform matrix from x_c in inertial frame to the origin of
        # inertial frame
        r_temp = [zeros(Float64, 3); -bs[i].x_c]
        r_temp = bs[i].Xb_to_i*r_temp
        r_temp = [zeros(Float64, 3); bs[i].x_i + r_temp[4:6]]
        Xic_to_i = TransMatrix(r_temp)
        # transform gravity force
        f_g = bs[i].Xb_to_i'*inv(Xic_to_i')*f_g
        # external force described in inertial coord
        f_ex = zeros(Float64, 6)
        f_ex = bs[i].Xb_to_i*f_ex
        # add up
        p_total[6i-5:6i] = p_temp - (f_g + f_ex)
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
        A_total[6i-5:6i, 6i-5:6i] = eye(Float64, 6)
        # fill in child joint blocks except for those body whose nchild=0
        for child_count = 1:bs[i].nchild
            chid = bs[i].chid[child_count]
            A_total[6i-5:6i, 6chid-5:6chid] = - bs[chid].Xp_to_b
        end
    end

    # collect all together
    return -p_total + A_total*sys.S_total*τ_total
end

#-------------------------------------------------------------------------------
function HERKFuncGT(bs::Vector{SingleBody}, sys::System)
"""
    HERKFuncGT constructs the input function GT for HERK method.
    It returns the force constraint matrix acting on Lagrange multipliers.
"""
    A_total = zeros(Float64, sys.ndof, sys.ndof)
    # construct A_total to take in parent-child hierarchy
    for i = 1:sys.nbody
        # fill in parent joint blocks
        A_total[6i-5:6i, 6i-5:6i] = eye(Float64, 6)
        # fill in child joint blocks except for those body whose nchild=0
        for child_count = 1:bs[i].nchild
            chid = bs[i].chid[child_count]
            A_total[6i-5:6i, 6chid-5:6chid] = - (bs[chid].Xp_to_b)'
        end
    end
    return A_total*sys.T_total
end

#-------------------------------------------------------------------------------
function HERKFuncG(bs::Vector{SingleBody}, sys::System)
"""
    HERKFuncG constructs the input function G for HERK method.
    It returns the motion constraint matrix acting on all body's velocity.
    These constraints arise from body velocity relation in each body's local
    body coord, for example if body 2 and 3 are connected then:
       v(3) = vJ(3) + X2_to_3*v(2)
"""
    B_total = zeros(Float64, sys.ndof, sys.ndof)
    # construct B_total to take in parent-child hierarchy
    for i = 1:sys.nbody
        # fill in child body blocks
        B_total[6i-5:6i, 6i-5:6i] = eye(Float64, 6)
        # fill in parent body blocks except for those body whose pid=0
        if bs[i].pid != 0
            pid = bs[i].pid
            B_total[6i-5:6i, 6pid-5:6pid] = - bs[i].Xp_to_b
        end
    end
    return (sys.T_total')*B_total
end

#-------------------------------------------------------------------------------
function HERKFuncgti(sys::System, t::T) where T <: AbstractFloat
"""
    HERKFuncgti returns all the collected prescribed active motion of joints
    at given time.
"""


end









end
