module TimeMarching

# export
export HERK!, Soln, BlockLU

# use registered packages
using DocStringExtensions

# import self-defined modules
using ..ConstructSystem
using ..SpatialAlgebra
using ..UpdateSystem
using ..Utils

#-------------------------------------------------------------------------------
function HERKScheme(name::String)
"""
HERKScheme provides a set of HERK coefficients, expressed in Butcher table form.
    A: Runge-Kutta matrix in Butcher tableau    c | A
    b: weight vector in Butcher tableau         ------
    c: node vector in Butcher tableau             | b
    st: stage number
    p: method order of accuracy
"""
    # Scheme A of 2-stage HERK in Liska's paper
    if name == "Liska"
        A = [0.0 0.0 0.0;
             0.5 0.0 0.0;
             √3/3 (3.0-√3)/3 0.0]
        b = [(3.0+√3)/6, -√3/3, (3.0+√3)/6]
        c = [0.0, 0.5, 1.0]
        st = 3
        p = 2
    # Brasey-Hairer 3-Stage HERK, table 2
    elseif name == "BH3"
        A = [0.0 0.0 0.0;
             1.0/3 0.0 0.0;
             -1.0 2.0 0.0]
        b = [0.0, 0.75, 0.25]
        c = [0.0, 1.0/3, 1.0]
        st = 3
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
        st = 5
        p = 4
    else
        error("This HERK scheme doesn't exist now.")
    end
    # modify for last stage
    A = [A; b']
    c = [c; 1.0]

    return A, b, c, st
end

#-------------------------------------------------------------------------------
function HERK!(sᵢₙ::Soln{T}, bs::Vector{SingleBody}, js::Vector{SingleJoint},
    sys::System) where T <: AbstractFloat
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
    @get sys.num_params (scheme, tol)

    # pick sheme parameters
    A, b, c, st = HERKScheme(scheme)
    if st != sys.num_params.st error("Scheme stage not correctly specified") end

    qJ_dim = sys.ndof
    λ_dim =sys.ncdof_HERK

    # pointer to pre-allocated array
    @get sys.pre_array (qJ, vJ, v, v̇, λ, v_temp, Mᵢ₋₁, fᵢ₋₁, GTᵢ₋₁, Gᵢ, gtiᵢ,
        lhs, rhs)

    # stage 1
    tᵢ₋₁ = sᵢₙ.t; tᵢ = sᵢₙ.t;; dt = sᵢₙ.dt
    qJ[1,:] = sᵢₙ.qJ
    v[1,:] = sᵢₙ.v
    # update vJ using v
    bs, js, sys, vJ[1,:] = UpdateVelocity!(bs, js, sys, v[1,:])

    # stage 2 to st+1
    for i = 2:st+1
        # time of i-1 and i
        tᵢ₋₁ = tᵢ
        tᵢ = sᵢₙ.t + dt*c[i]
        # initialize qJ[i,:]
        qJ[i,:] = sᵢₙ.qJ
        # calculate M, f and GT at tᵢ₋₁
        Mᵢ₋₁ = HERKFuncM(sys)
        fᵢ₋₁ = HERKFuncf(bs, js, sys)
        GTᵢ₋₁ = HERKFuncGT(bs, sys)
        # advance qJ[i,:]
        for k = 1:i-1
            qJ[i,:] += dt*A[i,k]*view(vJ,k,:)
        end
        # use new qJ to update system position
        bs, js, sys = UpdatePosition!(bs, js, sys, qJ[i,:])
        # calculate G and gti at tᵢ
        Gᵢ = HERKFuncG(bs, sys)
        gtiᵢ = HERKFuncgti(js, sys, tᵢ)
        # construct lhs matrix
        lhs = [ Mᵢ₋₁ GTᵢ₋₁; Gᵢ zeros(T,λ_dim,λ_dim) ]
        # the accumulated v term on the right hand side
        v_temp = sᵢₙ.v
        for k = 1:i-2
            v_temp += dt*A[i,k]*view(v̇,k,:)
        end
        # construct rhs
        rhs = [ fᵢ₋₁; -1./(dt*A[i,i-1])*(Gᵢ*v_temp + gtiᵢ) ]
######### use Julia's built in "\" operator for now
        # solve the eq
        x = lhs \ rhs
        # x = BlockLU(lhs, rhs, qJ_dim, λ_dim)
        # apply the solution
        v̇[i-1,:] = x[1:qJ_dim]
        λ[i-1,:] = x[qJ_dim+1:end]
        # advance v[i,:]
        v[i,:] = sᵢₙ.v
        for k = 1:i-1
            v[i,:] += dt*A[i,k]*view(v̇,k,:)
        end
        # update vJ using updated v
        bs, js, sys, vJ[i,:] = UpdateVelocity!(bs, js, sys, v[i,:])
# println("v = ", v[i,:])
    end

    # use norm(v[st+1,:]-v[st,:]) to determine next timestep
    sₒᵤₜ = Soln(tᵢ) # init struct
    sₒᵤₜ.dt = sᵢₙ.dt*(tol/norm(view(v,st+1,:)-view(v,st,:)))^(1/3)
    sₒᵤₜ.t = sᵢₙ.t + sᵢₙ.dt
    sₒᵤₜ.qJ = view(qJ, st+1, :)
    sₒᵤₜ.v = view(v, st+1, :)
    sₒᵤₜ.v̇ = view(v̇, st, :)
    sₒᵤₜ.λ = view(λ, st, :)

    return  sₒᵤₜ, bs, js, sys
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
    # pointer to pre-allocated array
    @get sys.pre_array (p_total, τ_total, p_bias, f_g, f_ex, r_temp,
        Xic_to_i, A_total)

    # compute bias force, gravity and external force
    for i = 1:sys.nbody
        # bias force
        p_bias = Mfcross(bs[i].v, (bs[i].inertia_b*bs[i].v))
        # gravity in inertial center coord
        f_g = bs[i].mass*[zeros(Float64, 3); sys.g]
        # get transform matrix from x_c in inertial frame to the origin of
        # inertial frame
        r_temp = [zeros(Float64, 3); -bs[i].x_c]
        r_temp = bs[i].Xb_to_i*r_temp
        r_temp = [zeros(Float64, 3); -bs[i].x_i + r_temp[4:6]]
        Xic_to_i = TransMatrix(r_temp)
        # transform gravity force
        f_g = bs[i].Xb_to_i'*inv(Xic_to_i')*f_g
        # external force described in inertial coord
        f_ex = zeros(Float64, 6)
        f_ex = bs[i].Xb_to_i*f_ex
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
    # pointer to pre-allocated array
    @get sys.pre_array (A_total,)

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
    # pointer to pre-allocated array
    @get sys.pre_array (B_total,)

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
function HERKFuncgti(js::Vector{SingleJoint}, sys::System, t::T) where
    T <: AbstractFloat
"""
    HERKFuncgti returns all the collected prescribed active velocity of joints
    at given time.
"""
    # pointer to pre-allocated array
    @get sys.pre_array (v_gti, va_gti)

    # give actual numbers from calling motion(t)
    for i = 1:sys.na
        jid = sys.kinmap[i,1]
        dofid = sys.kinmap[i,2]
        _, va_gti[i] = js[jid].joint_dof[dofid].motion(t)
    end

    v_gti[sys.udof_a] = va_gti
    return -(sys.T_total')*v_gti
end

#-------------------------------------------------------------------------------
function BlockLU(H::Array{T,2}, b::Vector{T}, qJ_dim::Int,
    λ_dim::Int) where T
"""
    BlockLU solve the system H*xy = b using Schur complement reduction
    [A  B₁ᵀ] * [x] = [f]
    [B₂  -C]   [y]   [g]
    -------   ---   ---
       H    * xy  =  b
    By computing the Schur complement S = -B₂*inv(A)*B₁ᵀ-C, the original system
    of equations is transferred to
    [A  B₁ᵀ] * [x] = [       f       ]
    [0   S ]   [y]   [g - B₂*inv(A)*f]
    Using lufact to solve the lower part y first, then get x by plugging in y.
"""
    # set pointers
    A = view(H, 1:qJ_dim, 1:qJ_dim)
    B₁ᵀ = view(H, 1:qJ_dim, qJ_dim+1:qJ_dim+λ_dim)
    B₂ = view(H, qJ_dim+1:qJ_dim+λ_dim, 1:qJ_dim)
    C = - view(H, qJ_dim+1:qJ_dim+λ_dim, qJ_dim+1:qJ_dim+λ_dim)
    f = view(b, 1:qJ_dim)
    g = view(b, qJ_dim+1:qJ_dim+λ_dim)
    # compute Schur complement S
    S = - B₂*inv(A)*B₁ᵀ - C
    # compute y first, then compute x by substitute in y
    y = S \ (g - B₂*inv(A)*f)
    return [A \ (f - B₁ᵀ*y); y]
end




end
