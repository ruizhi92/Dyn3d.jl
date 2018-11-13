@testset "FluidInteraction" begin
    @testset "DetermineNP" begin
        @test DetermineNP(1, 0.02) == 201
        @test DetermineNP(1, 0.01) == 401
        @test DetermineNP(2, 0.02) == 101
        @test DetermineNP(3, 0.02) == 65
        @test DetermineNP(4, 0.02) == 49
    end

    # Construct a BodyDyn structure
    include(Pkg.dir("Dyn3d")*"/test/config_body.jl")
    bodys, joints, system = BuildChain(config_bodys, config_joints, config_system)
    bd = BodyDyn(bodys, joints, system)
    bd, soln = InitSystem!(bd)
    herkbody = Dyn3d.HERKBody(system.num_params,HERKFuncM, HERKFuncGT, HERKFuncG,
                (HERKFuncf,HERKFuncgti), (UpdatePosition!,UpdateVelocity!))

    # Construct fluid
    include(Pkg.dir("Dyn3d")*"/test/config_fluid.jl")

    # create fluid-body interface
    bgs = GenerateBodyGrid(bd; np=DetermineNP(nbody, Δx))
    @test bgs[1].np == 201
    @test length(bgs[1].points) == 201

    # cut out in 2d
    bgs = CutOut2d(bd,bgs)
    @test bgs[1].np == 51
    @test length(bgs[1].points) == 51

    # check intial body points position
    bgs = AcquireBodyGridKinematics(bd,bgs)
    @test isapprox(hcat(bgs[1].q_i...)'[:,1],1:0.02:2)

    # advance body solver for one step
    soln.dt = 0.01
    soln, bds = herkbody(soln, bd; _isfixedstep=true, _outputmode=true)
    NS = 3

    # check body points position & velocity
    bkins = Vector{Array{Float64,2}}(NS)
    for k = 1:NS
        bgs = AcquireBodyGridKinematics(bds[k],bgs)
        coord = hcat(bgs[1].q_i...)'[:,[1,2]]
        motion = hcat(bgs[1].v_i...)'[:,[1,2]]
        for i = 2:length(bgs)
            coord = [coord[1:end-1,:]; hcat(bgs[i].q_i...)'[:,[1,2]]]
            motion = [motion[1:end-1,:]; hcat(bgs[i].v_i...)'[:,[1,2]]]
        end
        bkins[k] = [coord motion]
    end
    @test all(bkins[3][:,2] .≈ 0.99995)
    @test all(bkins[3][:,4] .≈ -0.01)

    # create fluid solver and advance for one step
    coord_init = bkins[1][:,1:2]
    X̃ = VectorData(coord_init)
    sys = Systems.NavierStokes((nx,ny),Re,Δx,Δt,U∞ = U∞, X̃ = X̃, isstore = true, isstatic = false)
    u = w₀
    f = VectorData(X̃)
    fs = [[f.u f.v] for i=1:NS]
    ifherk_sc2d = Whirl.IFHERK_sc2d(u,f,sys.Δt,
                    (t,u) -> Systems.plan_intfact(t,u,sys),
                    (u,t,coord) -> Whirl.plan_constraints(u,t,sys,coord),
                    ((u,t) -> Whirl.r₁(u,t,sys),
                     (u,t,motion) -> Whirl.r₂(u,t,sys,motion)),
                    coord_init,
                    tol=1e-3,rk=Whirl.TimeMarching.RK31,
                    isstored=true,isstaticconstraints=false)
    t, u, f, fs = ifherk_sc2d(t,u,bkins)

    # integrate body force
    for i = 1:size(fs[1],1)
        bgs[1].f_ex3d[i][[1,2]] = fs[3][i,:]*Δx^2
    end
    bgs = IntegrateBodyGridDynamics(bd,bgs)

    # there's no rotation for this special case
    @test isapprox(sum(fs[3][:,1])*Δx^2, bgs[1].f_ex6d[4]; atol=1e-7)
    @test isapprox(sum(fs[3][:,2])*Δx^2, bgs[1].f_ex6d[5]; atol=1e-7)

end
