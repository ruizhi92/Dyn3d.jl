var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Dyn3d.jl-1",
    "page": "Home",
    "title": "Dyn3d.jl",
    "category": "section",
    "text": "A 2d/3d rigid body dynamics solverThe main goal of this repository is to construct a rigid body-joint system and solve forward/backward dynamics problem on it. This package is functioned through:constructing 2d polygon shape rigid bodies and allow motion in 2d/3d space\nconnecting bodies by joints which has 6 degree of freedoms for each\nsolving motions on unconstrained(passive) degrees of freedom of joints\nsolving forces on constrained degrees of freedom, allowing active motion\nplotting/making gif in Julia or making movies in MatlabTo solve a rigid body dynamics problem, this package express the dynamics using 6D spatial vector developed by Roy Featherstone[1]. The governing equations are formulated to fit in half explicit Runge-Kutta method on index-2 differential equations[2]. Constrained forces on joints are represented in Lagrange multiplier terms and solved together with motions of all degrees of freedom.Based on the calculation of dynamical systems, Dyn3d.jl is also used to simulate fluid-structure interaction(FSI) problems together with package Whirl.jl for strongly coupled method. Notebook example is provided in notebook folder. Fully coupled method taking advantage of both Dyn3d.jl and Whirl.jl is implemented in package FSI.jl.(Image: )"
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "This package currently requires Julia 0.6. To install, simply download this Github repository, find the location of this repository expressed in Julia byjulia> Dyn3d_dir = Pkg.dir(\"Dyn3d\")and then setup a symbolic link in shell followingshell$ sudo ln -s actual_address Dyn3d_dirThe plots in this documentation are generated using Plots.jl. You might want to install that too to follow the examples.If you have trouble in setting up the symbolic like to directory of Dyn3d.jl, a simple alternative solution is:include(path*\"Dyn3d.jl\")"
},

{
    "location": "index.html#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "[1]: Featherstone, Roy. Rigid body dynamics algorithms. Springer, 2014.[2]: Brasey, Valérie, and Ernst Hairer. \"Half-explicit Runge–Kutta methods for differential-algebraic systems of index 2.\" SIAM Journal on Numerical Analysis 30, no. 2 (1993): 538-552."
},

{
    "location": "manual/construct_system.html#",
    "page": "Construct body-joint system",
    "title": "Construct body-joint system",
    "category": "page",
    "text": ""
},

{
    "location": "manual/construct_system.html#Construct-body-joint-system-1",
    "page": "Construct body-joint system",
    "title": "Construct body-joint system",
    "category": "section",
    "text": "Constructing the body-joint system requires filling in a lot of detailed information, like body-joint hierarchy information, body inertia, joint degree of freedoms type etc. In this section we introduce the construction of the system in two aspects:User-oriented set up(front end)\nCode-oriented construction(back end)"
},

{
    "location": "manual/construct_system.html#User-oriented-set-up-1",
    "page": "Construct body-joint system",
    "title": "User-oriented set up",
    "category": "section",
    "text": "To set up a body-joint system, one needs to know how to set body properties and joint properties, also the hierarchy information. In this section we introduce how to set up the body-joint system, together with other information needed to enclosure the problem we\'re interested in.The set up of the problem generally has 3 parts: set up body, set up joint, and set up system information like the dimension of this problem, gravity, numerical parameters etc. We will use these information provided by the user to construct every single body, every single joint, and finally connect them with hierarchy information. Some examples of providing user-set-up information are listed in src/config_files for user convenience. Here we take the 2dFall.jl as an example to show the procedure of setting up the problem.ConfigDataType"
},

{
    "location": "manual/construct_system.html#Code-oriented-construction-1",
    "page": "Construct body-joint system",
    "title": "Code-oriented construction",
    "category": "section",
    "text": "Before running any dynamics, we need to construct the body-joint system, connecting bodies by joints in the correct hierarchy using supplied configuration information. In order to do that, three key functions are used:AddBody\nAddJoint\nAssembleSystemwhich is responsible for add a body to the body system, add a joint to the joint system and assemble the body-joint system and fill in extra hierarchy information. These three functions are also bundled in the correct sequential order inBuildChain(cbs, cjs, csys)In order to solve for the rigid body-joint chain system, we choose body velocity v and joint displacement qJ as main variables in the time marching scheme. We also need to construct a structure that stores and updates the body-joint chain\'s intermediate variables such as all kinds of transformation matrices X, body position in the inertial space x_i and so on. All information of the body-joint system is bundled into a BodyDyn structure."
},

{
    "location": "manual/construct_system.html#Dyn3d.ConfigDataType.ConfigBody",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConfigDataType.ConfigBody",
    "category": "type",
    "text": "ConfigBody(nbody::Int,nverts::Int,verts::Matrix{Float64},ρ::Float64)\n\nSet up configuration information for the a single body in the body system. Here we assume that all bodies has the same shape if more than one body exists. A single body is an infinitely thin body in y direction. It must have polygon shape and is described in z-x space. For example if we describe a rectangle in z-x space, for 3d problem it\'s just fine. For 2d problem in x-y space, this rectangle has a projection  of a line in x-y space. The vertices local coordinates are described in clockwise direction as a convention.\n\nFields\n\nnbody: Number of bodies in total\nnverts: Number of vertices for one body\nverts: Polygon vertices coordinates starting from the left bottom vert and   going in clockwise direction. Each line describes the (z,x) coordinate of   this vertice. Usually the z-dimesion has unit length 1 for 2-d problem.\nρ: Density of this body in mass per area\n\nBody setup\n\n     ^(y)\n     |\n     |----->(x)\n    /\n (-z)\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConfigDataType.ConfigJoint",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConfigDataType.ConfigJoint",
    "category": "type",
    "text": "ConfigJoint(njoint,joint_type,shape1,shape2,body1,joint_dof,qJ_init)\n\nSet up configuration information for a single joint. A joint allows one/multiple degree of freedoms from 1 to 6.\n\nFields\n\nnjoint: Int, total number of joints for this body-joint system\njoint_type: String, allows \"revolute\", \"prismatic\", \"cylindrical\", \"planar\",   \"spherical\", \"free\" and \"custom\". Detailed information is described in JointType   section\nshape1: Vector{Float64} with 6 elements. It describes the location of this   joint in its parent body coordinate. If shape1 is used on the first joint,   then it\'s the orientation of this first joint to inertial system. It is written   with [θx, θy, θz, x, y, z].\nshape2: Vector{Float64} with 6 elements. It describes the location of this   joint in its child body coordinate. Normally, shape2 is zeros(Float64,6) if   there\'s no distance gap or angle gap between two adjacent bodies.\nbody1: the parent body id of this joint\njoint_dof: Vector{Float64}. The size of joint_dof depends on the number of   specified dof, refers to type Dof.\nqJ_init: Vector{Float64}. It is the initial angle/displacement of this joint   with respect to its parent body. The size should be the same as the number of   dof specified.\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConfigDataType.ConfigSystem",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConfigDataType.ConfigSystem",
    "category": "type",
    "text": "ConfigSystem(ndim::Int,gravity::Vector{Float64},num_params::NumParams)\n\nAdditional system information to define this problem.\n\nFields\n\nndim: Dimension of this problem. Choices are 2 or 3\ngravity: Non-dimensional gravity. It is in [x,y,z] direction. So if we\'re   describing a 2d problem with gravity pointing downward, it should be [0.,-1.,0.]\nnum_params: Refer to type NumParams\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConfigDataType.Dof",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConfigDataType.Dof",
    "category": "type",
    "text": "Dof(dof_id::Int,dof_type::String,stiff::Float64,damp::Float64,motion::Motions)\n\nSet up a single degree of freedom(dof) information in a joint. A joint may have a maximum 6 dofs and minimum 1 dof(either active or passive). Here we don\'t allow it to have 0 since there\'s no reason to do this. If we want the parent body and child body to have no relative motion(i.e. they\'re rigidly connected together), we can set the second joint has only one dof and this dof has active \"hold\" motion.\n\nFields\n\ndof_id: Choose from 1 to 6 and is corresponding to [Ox, Oy, Oz, x, y, z].\ndof_type: \"passive\" or \"active\"\nstiff: Non-dimensional stiffness of a string associated with this degree of   freedom. Only has effect on solving the system when this dof is passive.\ndamp: Similar to stiff, this assigns the damping coefficient of this dof.\nmotion: Defines the motion of this dof. Refer to type Motion\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConfigDataType.Motions",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConfigDataType.Motions",
    "category": "type",
    "text": "Motions(type::String, parameters::Vector{Float64})\n\nA structure representing active motion of a joint, allowing different types of motion, can be time-dependent.\n\nFields\n\nmotion_type: Allow choices of \"hold\", \"velocity\", \"oscillatory\", \"ramp_1\", \"ramp_2\"\nmotion_params: Numerical parameters provided to describe the motion\n\nConstructors\n\nMotions(): Provide no active motion\nMotions(\"hold\", [qJ]): Hold the joint at constant qJ and vJ=0\nMotions(\"velocity\",[qJ,vJ]): Specify constant vJ of this joint with initial angle qJ\nMotions(\"oscillatory\",[amp,freq,phase]) specify a oscillatory motion through   qJ = amp*cos(2*freq*t+phase)\nMotions(\"ramp_1\",[a,t₁,t₂,t₃,t₄]): Describes a ramp motion in [1]\nMotions(\"ramp_2\",[a]): Describes a decelerating motion with an initial velocity\n\n[1]: Eldredge, Jeff, Chengjie Wang, and Michael Ol. \"A computational study of a canonical pitch-up,\n\npitch-down wing maneuver.\" In 39th AIAA fluid dynamics conference, p. 3687. 2009.\n\nFunctions\n\n(m::Motions)(t): Evalutes the motion type at time t, return joint angle qJ and   velocity vJ\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConfigDataType.NumParams",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConfigDataType.NumParams",
    "category": "type",
    "text": "NumParams(tf::Float64,dt::Float64,scheme::String,st::Int,tol::Float64)\n\nNumerical parameters needed for the time marching scheme.\n\nFields\n\ntf: The end time of this run\ndt: Time step size\nscheme: Applies the implicit Runge-kutta method of different coefficient, choices   are \"Liska\"(2nd order), \"BH3\"(3rd order), \"BH5\"(4th order), \"Euler\"(1st order),   \"RK2\"(2nd order), \"RK22\"(2nd order).\nst: The number of stages of this RK scheme\ntol: Tolerance used in time marching for adptive time step\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConstructSystem.Soln",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConstructSystem.Soln",
    "category": "type",
    "text": "This module construct the body-joint system by:     1. AddBody     2. AddJoint     3. AssembleSystem\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConstructSystem.AddBody-Tuple{Int64,Dyn3d.ConfigDataType.ConfigBody}",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConstructSystem.AddBody",
    "category": "method",
    "text": "AddBody(id::Int, cf::ConfigBody)\n\nUsing ConfigBody information, add one body to the body system bs\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConstructSystem.AddJoint-Tuple{Int64,Dyn3d.ConfigDataType.ConfigJoint}",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConstructSystem.AddJoint",
    "category": "method",
    "text": "AddJoint(id::Int, cf::ConfigJoint)\n\nUsing ConfigJoint information, add one joint to the joint system js\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConstructSystem.AssembleSystem-Tuple{Array{Dyn3d.ConstructSystem.SingleBody,1},Array{Dyn3d.ConstructSystem.SingleJoint,1},Dyn3d.ConstructSystem.System}",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConstructSystem.AssembleSystem",
    "category": "method",
    "text": "AssembleSystem(bs::Vector{SingleBody}, js::Vector{SingleJoint}, sys::System)\n\nWith body system bs and joint system js, fill in extra hierarchy information about the whole system, connects bodies and joints.\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Dyn3d.ConstructSystem.BuildChain-Tuple{Array{Dyn3d.ConfigDataType.ConfigBody,1},Array{Dyn3d.ConfigDataType.ConfigJoint,1},Dyn3d.ConfigDataType.ConfigSystem}",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConstructSystem.BuildChain",
    "category": "method",
    "text": "BuildChain(cbs::Vector{ConfigBody}, cjs::Vector{ConfigJoint}, csys::ConfigSystem)\n\nPut AddBody, AddJoint and AssembleSystem in sequential order in a single function\n\n\n\n"
},

{
    "location": "manual/construct_system.html#Methods-1",
    "page": "Construct body-joint system",
    "title": "Methods",
    "category": "section",
    "text": "Modules = [ConfigDataType, ConstructSystem]\nOrder   = [:type, :function]"
},

{
    "location": "manual/construct_system.html#Index-1",
    "page": "Construct body-joint system",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"construct_system.md\"]"
},

{
    "location": "manual/fluid_interaction.html#",
    "page": "Fluid-Structure Interaction",
    "title": "Fluid-Structure Interaction",
    "category": "page",
    "text": ""
},

{
    "location": "manual/fluid_interaction.html#Fluid-Structure-Interaction-1",
    "page": "Fluid-Structure Interaction",
    "title": "Fluid-Structure Interaction",
    "category": "section",
    "text": "DocTestSetup = quote\nusing Whirl\nusing Dyn3d\nusing Plots\nenddefddt1fracmathrmd1mathrmdta = 1Dyn3d can be used to solve rigid body dynamics, to simulate fluid-structure interaction, we use Whirl.jl branch stronglycoupled2d to solve for fluid dynamics. To make fluid and body dynamics coupling with each other, we also need to construct a series of functions to act as interface between them, and choose appropriate coupling schemes.Strongly coupled method is finished, and fully coupled method is under construction. Example of working fluid-structure interaction notebook is also provided under /notebook."
},

{
    "location": "manual/fluid_interaction.html#Strongly-coupled-method-1",
    "page": "Fluid-Structure Interaction",
    "title": "Strongly coupled method",
    "category": "section",
    "text": "Coupling scheme is refered to Wang and Eldredge JCP, detailes are skipped here. The main idea is to use fluid forces on body as the iterating variable, run the body solver Dyn3d with proposed fluid force as external force for the body-joint system, get updated body coordinates and velocities. Then use these information as boundary condition for fluid solver Whirl, get updated fluid force on body from Lagrange multipliers. If the proposed force and updated force are close enough, this timestep is said to be converged. Otherwise we use a relaxation scheme to calculate a new proposed force, iterates until it converge.There are several things that needs to be carefully dealt with, because we are using two packages together and variables need to be consistent from one to another."
},

{
    "location": "manual/fluid_interaction.html#Scaling-1",
    "page": "Fluid-Structure Interaction",
    "title": "Scaling",
    "category": "section",
    "text": "For example gravity in body solver should be non-dimensionalized to 00 -10 00 instead of 00 -98 00. Scaling is done through:mass ratio m^* = fracrho_bLrho_f = fracrho_brho_f LReynolds number Re = fracrho_f U_infty Lmutorsion stiffness k^* = frackrho_f U_infty^2 L^2bending stiffness c^* = fraccrho_f U_infty L^3gravity g^* = fracg LU_infty^2 if uniform flow is non-zero"
},

{
    "location": "manual/fluid_interaction.html#Fully-coupled-method-1",
    "page": "Fluid-Structure Interaction",
    "title": "Fully coupled method",
    "category": "section",
    "text": "Under construction"
},

{
    "location": "manual/fluid_interaction.html#Dyn3d.FluidInteraction.BodyGrid",
    "page": "Fluid-Structure Interaction",
    "title": "Dyn3d.FluidInteraction.BodyGrid",
    "category": "type",
    "text": "BodyGrid(bid::Int, np::Int, points, q_i, f_ex3d, f_ex6d)\n\nDesign this structure to contain body points coord, motion and forces used in fluid-structure interaction\n\nFields\n\nbid: body id in the joint-body chain\nnp: number of grid points on this body\npoints: (x,y,z) coordinates of all body points in local body frame\nq_i: (x,y,z) coordinates of all body points in inertial frame\nv_i: velocity of all body points in inertial frame\nf_ex3d: external force(fluid force) on all body points in inertial frame\nf_ex6d: f_ex3d integrated through all body points on one body and described in 6d spatial vector form\n\nConstructors\n\nBodyGrid(bid,np,points): initialize q_i, v_i, f_ex3d to be Vector of zeros(3),                            f_ex6d to be zeros(6)\n\n\n\n"
},

{
    "location": "manual/fluid_interaction.html#Dyn3d.FluidInteraction.AcquireBodyGridKinematics-Tuple{Dyn3d.ConstructSystem.BodyDyn,Array{Dyn3d.FluidInteraction.BodyGrid,1}}",
    "page": "Fluid-Structure Interaction",
    "title": "Dyn3d.FluidInteraction.AcquireBodyGridKinematics",
    "category": "method",
    "text": "AcquireBodyGridKinematics(bd::BodyDyn, bgs::Vector{BodyGrid})\n\nGiven updated bd structure, which contains 3d bs[i].x_i in inertial frame and 6d bs[i].v of each body in the body local frame, return 3d linear q_i and v_i of each body point in the inertial frame.\n\n\n\n"
},

{
    "location": "manual/fluid_interaction.html#Dyn3d.FluidInteraction.CutOut2d-Tuple{Dyn3d.ConstructSystem.BodyDyn,Array{Dyn3d.FluidInteraction.BodyGrid,1}}",
    "page": "Fluid-Structure Interaction",
    "title": "Dyn3d.FluidInteraction.CutOut2d",
    "category": "method",
    "text": "CutOut2d(bd::BodyDyn,bgs::Vector{BodyGrid})\n\nThis function need to be called only once after GenerateBodyGrid for 2d case of flat plates.\n\nIn Dyn3d, bodies are constructed by quadrilateral/triangles(not lines) in z-x plane for both 2d/3d cases. In Whirl, fluid in 2d cases are constructed in x-y plane. Thus to describe plates as lines in x-y space, we cut out the info on the other sides of the plate. Note that verts are formulated in clockwise direction, with the left-bottom corner as origin.\n\n\n\n"
},

{
    "location": "manual/fluid_interaction.html#Dyn3d.FluidInteraction.DetermineNP-Tuple{Int64,Float64}",
    "page": "Fluid-Structure Interaction",
    "title": "Dyn3d.FluidInteraction.DetermineNP",
    "category": "method",
    "text": "DetermineNP(nbody::Int, Δx)\n\nRun this function before running GenerateBodyGrid, to determine number of points on a 2d body, in order to satisfy the desired number of points on the 1d body.\n\nnp = (# of points on 1d plate - 1)*4+1. So np=201 has 51 points(1 body), np=101 has 26 points(2 body), np=49 has 13 points(4 body), np=25 has 7 points(8 body), etc.\n\n\n\n"
},

{
    "location": "manual/fluid_interaction.html#Dyn3d.FluidInteraction.GenerateBodyGrid-Tuple{Dyn3d.ConstructSystem.BodyDyn}",
    "page": "Fluid-Structure Interaction",
    "title": "Dyn3d.FluidInteraction.GenerateBodyGrid",
    "category": "method",
    "text": "GenerateBodyGrid(bd::BodyDyn; np=101)\n\nGiven BodyDyn structure, where each body only consists of several verts(usually 4 for quadrilateral and 3 for triangle), return the verts position in inertial frame of given number of points np by interpolation, of all bodies in the system.\n\n\n\n"
},

{
    "location": "manual/fluid_interaction.html#Dyn3d.FluidInteraction.IntegrateBodyGridDynamics-Tuple{Dyn3d.ConstructSystem.BodyDyn,Array{Dyn3d.FluidInteraction.BodyGrid,1}}",
    "page": "Fluid-Structure Interaction",
    "title": "Dyn3d.FluidInteraction.IntegrateBodyGridDynamics",
    "category": "method",
    "text": "IntegrateBodyGridDynamics(bd::BodyDyn, bgs::Vector{BodyGrid})\n\nGiven external 3d linear fluid force f_ex of each body point contained in updated bgs structure, do intergral to return integrated 6d body force([torque,force]) exerting on the beginning of current body, desribed in inertial frame.\n\n\n\n"
},

{
    "location": "manual/fluid_interaction.html#Methods-1",
    "page": "Fluid-Structure Interaction",
    "title": "Methods",
    "category": "section",
    "text": "Modules = [FluidInteraction]\nOrder   = [:type, :function]"
},

{
    "location": "manual/fluid_interaction.html#Index-1",
    "page": "Fluid-Structure Interaction",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"fluid_interaction.md\"]"
},

]}
