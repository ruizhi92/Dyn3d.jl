@testset "JointType" begin
    kind_group = ["revolute", "prismatic", "cylindrical", "planar", "spherical"]
    for i = 1:length(kind_group)
        j = Dyn3d.ChooseJoint(kind_group[i])
        @test j.nudof == length(j.udof)
        @test j.ncdof == length(j.cdof)
        @test (6,j.nudof) == size(j.S)
        @test (6,j.ncdof) == size(j.T)
    end

end

# tests through every config files in Config_files folder
@testset "ConstructSystem" begin
    path = Pkg.dir("Dyn3d")*"/src/Config_files"
    names = readdir(path)

    for k = 1:length(names)
        if names[k][1] == '2' || names[k][1] == '3'
            # load config files
            include(path*'/'*names[k])

            # check ConfigDataType
            @test njoint == nbody
            if isdefined(:config_bodys)
                for i = 1:length(config_bodys)
                    @test config_bodys[i].nverts == size(config_bodys[i].verts,1)
                    @test length(config_joints[i].joint_dof) ==
                          length(config_joints[i].qJ_init)
                end
            else
                @test config_body.nverts == size(config_body.verts,1)
            end

            # check AddBody
            bodys = Vector{Dyn3d.SingleBody}(nbody)
            for i = 1:nbody
                bodys[i] = Dyn3d.AddBody(i, config_body)
            end

            # check AddJoint
            joints = Vector{Dyn3d.SingleJoint}(njoint)
            for i = 1:njoint
                joints[i] = Dyn3d.AddJoint(i, config_joints[i])
            end

            # check AssembleSystem
            system = Dyn3d.System(ndim, nbody, njoint, gravity, num_params)
            bodys, joints, system = Dyn3d.AssembleSystem!(bodys, joints, system)
        end
    end


end
