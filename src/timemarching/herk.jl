 """
    HERK(u,f,Δt,plan_intfact,B₁ᵀ,B₂,r₁,r₂;[tol=1e-3],[issymmetric=false],[rk::RKParams=RK31])
    HERKBody

    HERK is a half-explicit Runge-Kutta solver based on the
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

# Arguments

- `u` : example of state vector data
- `f` : example of constraint force vector data
- `Δt` : time-step size
- `plan_intfact` : constructor to set up integrating factor operator for `A` that
              will act on type `u` (by left multiplication) and return same type as `u`
- `plan_constraints` : constructor to set up the
- `B₁ᵀ` : operator acting on type `f` and returning type `u`
- `B₂` : operator acting on type `u` and returning type `f`
- `r₁` : operator acting on type `u` and `t` and returning `u`
- `r₂` : operator acting on type `u` and `t` and returning type `f`
"""

mutable struct HERKBody{FA,FB1,FB2,FR1,FR2,FP,FV}

    # numerical parameters
    rk :: RKParams
    tol :: Float64

    # left hand side matrix components
    A :: FA # operates on TU and returns TF
    B₁ᵀ :: FB1  # operates on TF and returns TU
    B₂ :: FB2   # operates on TU and returns TF
    r₁ :: FR1  # function of u and t, returns TU
    r₂ :: FR2  # function of t, returns TF

    # functions to update lhs matrix components
    UpP :: FP # update joint.qJ
    UpV :: FV # update body.v

    # # Saddle-point systems
    # S :: Vector{SaddleSystem}  # -B₂AB₁ᵀ
end

#-------------------------------------------------------------------------------
function (::Type{HERKBody})(num_params::NumParams, A::FA, B₁ᵀ::FB1, B₂::FB2,
                            rhs::Tuple{FR1,FR2}, up::Tuple{FP,FV},
                            ) where {FA,FB1,FB2,FR1,FR2,FP,FV}

    @get num_params (scheme, tol)
    rk = RKParams(scheme)

    return HERKBody{FA,FB1,FB2,FR1,FR2,FP,FV}(rk, tol, A, B₁ᵀ,
                B₂, rhs[1], rhs[2], up[1], up[2])
end

#-------------------------------------------------------------------------------
function Base.show(io::IO, scheme::HERKBody{FA,FB1,FB2,FR1,FR2,FP,FV}) where {FA,FB1,FB2,FR1,FR2,FP,FV}
    println(io, "Order-$(scheme.rk.st) HERK integrator.")
end

#-------------------------------------------------------------------------------
function (scheme::HERKBody{FA,FB1,FB2,FR1,FR2,FP,FV})(sᵢₙ::Soln{T}, bd::BodyDyn) where
        {T<:AbstractFloat,FA,FB1,FB2,FR1,FR2,FP,FV}

    @get scheme (rk, tol, A, B₁ᵀ, B₂, r₁, r₂, UpP, UpV)
    @get bd (bs, js, sys)
    @get rk (st, c, a)

    if st != sys.num_params.st error("Scheme stage not correctly specified") end

    bs, js, sys = bd.bs, bd.js, bd.sys

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
    bs, js, sys, vJ[1,:] = UpV(bs, js, sys, v[1,:])

    # stage 2 to st+1
    for i = 2:st+1
        # time of i-1 and i
        tᵢ₋₁ = tᵢ
        tᵢ = sᵢₙ.t + dt*c[i]
        # initialize qJ[i,:]
        qJ[i,:] = sᵢₙ.qJ
        # calculate M, f and GT at tᵢ₋₁
        Mᵢ₋₁ = A(sys)
        fᵢ₋₁ = r₁(bs, js, sys)
        GTᵢ₋₁ = B₁ᵀ(bs, sys)
        # advance qJ[i,:]
        for k = 1:i-1
            qJ[i,:] += dt*a[i,k]*view(vJ,k,:)
        end
        # use new qJ to update system position
        bs, js, sys = UpP(bs, js, sys, qJ[i,:])
        # calculate G and gti at tᵢ
        Gᵢ = B₂(bs, sys)
        gtiᵢ = r₂(js, sys, tᵢ)
        # construct lhs matrix
        lhs = [ Mᵢ₋₁ GTᵢ₋₁; Gᵢ zeros(T,λ_dim,λ_dim) ]
        # the accumulated v term on the right hand side
        v_temp = sᵢₙ.v
        for k = 1:i-2
            v_temp += dt*a[i,k]*view(v̇,k,:)
        end
        # construct rhs
        rhs = [ fᵢ₋₁; -1./(dt*a[i,i-1])*(Gᵢ*v_temp + gtiᵢ) ]
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
            v[i,:] += dt*a[i,k]*view(v̇,k,:)
        end
        # update vJ using updated v
        bs, js, sys, vJ[i,:] = UpV(bs, js, sys, v[i,:])
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
