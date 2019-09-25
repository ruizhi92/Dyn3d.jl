@testset "HERK" begin
    # preparation for testing HERK
    methods = ["Liska", "BH3", "BH5"]
    sts = [3, 3, 5]
    path = "../src/config_files/2dLink.jl"
    include(path)

    # begin time marching
    for k in 1:length(methods)
        # build body chain
        bodys, joints, system = BuildChain(config_bodys, config_joints, config_system)
        system.num_params.tf = 1e-3

        # use different methods
        system.num_params.scheme = methods[k]
        system.num_params.st = sts[k]

        # init system
        bd = BodyDyn(bodys, joints, system)
        bd, soln = InitSystem!(bd)

        herk = HERKBody(system.num_params,HERKFuncM, HERKFuncGT, HERKFuncG,
                        (HERKFuncf,HERKFuncgti), (UpdatePosition!,UpdateVelocity!))

        while soln.t < system.num_params.tf
            soln, bd = herk(soln, bd, _isfixedstep=true)
        end
    end
end

# @testset "BlockLU" begin
#     n = 30
#     t = rand(1:n)
#     A = rand(n,n)
#     b = rand(n)
#     @test BlockLU(A, b, t, n-t) ≈ A\b
# end

@testset "2d plate falling with initial velocity" begin
    # problem dimension
    ndim = 2
    # numerical params
    tf = 2
    dt = 1e-3
    scheme = "Liska"
    st = 3
    tol = 1e-4
    num_params = NumParams(tf, dt, scheme, st, tol)
    # gravity
    gravity = [0., -1.0, 0.]
    # set up system config info
    config_system = ConfigSystem(ndim, gravity, num_params)
    # set up bodys
    nbody = 1
    config_body = ConfigBody(nbody, 4,
       [0. 0.; 1. 0.; 1. 1.0/nbody; 0. 1.0/nbody], 1.0)
    config_bodys = fill(config_body, nbody)
    # set up joints
    njoint = nbody
    config_joints = Vector{ConfigJoint}(undef,njoint)
    # set the first passive joint with no stiff and damp
    dof_1 = Dof(5, "passive", 0., 0., Motions())
    # set initial downward velocity
    v₀ = 0.5
    config_joints[1] = ConfigJoint(njoint, "custom_prismatic_in_y",
        [0.,0.,0.,0.5,1.5,0.], zeros(Float64,6), 0, [dof_1], [0.], [v₀])

    # Build joint-body chain
    bodys, joints, system = BuildChain(config_bodys, config_joints,
                                       config_system)
    bd = BodyDyn(bodys, joints, system)

    # Initialize system state
    bd, soln = InitSystem!(bd)

    # init soln structure
    solns = (Soln)[]
    push!(solns, soln)

    # init VertsHistory struct
    vs = []
    push!(vs, VertsHistory(system.nbody, bd.bs));

    # Set up HERKBody object
    herk = HERKBody(system.num_params,HERKFuncM, HERKFuncGT, HERKFuncG,
                    (HERKFuncf,HERKFuncgti), (UpdatePosition!,UpdateVelocity!))

    # Time Marching
    idx = 0
    while soln.t < tf
        # advance one step
        soln, bd = herk(soln, bd, _isfixedstep=true)

        # record soln and verts_i info
        push!(solns, soln)
        push!(vs, VertsHistory(system.nbody, bd.bs))

        # print progress
        idx += 1
    end

    # check end position
    @test isapprox(-0.5*(tf+dt)^2 + v₀*(2+dt),solns[end].qJ[5],atol=1e-10)
    # check end velocity
    @test isapprox(v₀-(tf+dt),solns[end].v[5],atol=1e-10)


end
