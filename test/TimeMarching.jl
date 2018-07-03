@testset "HERK" begin
    # preparation for testing HERK
    methods = ["Liska", "BH3", "BH5"]
    sts = [3, 3, 5]
    path = Pkg.dir("Dyn3d")*"/src/config_files/2dLink.jl"
    include(path)
    bodys, joints, system = BuildChain(config_bodys, config_joints, config_system)
    system.num_params.tf = 1e-3
    # init system
    bodys, joints, system, soln = InitSystem!(bodys, joints, system)

    # begin time marching
    for k in 1:length(methods)
        system.num_params.scheme = methods[k]
        system.num_params.st = sts[k]
        while soln.t < system.num_params.tf
            soln, bodys, joints, system = HERK!(soln, bodys, joints, system)
        end
    end
end

@testset "BlockLU" begin
    n = 30
    t = rand(1:n)
    A = rand(n,n)
    b = rand(n)
    @test BlockLU(A, b, t, n-t) ≈ A\b
end
