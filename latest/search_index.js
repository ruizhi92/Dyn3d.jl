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
    "text": "A 2d/3d rigid body dynamics solverThe main goal of this repository is to construct a rigid body-joint system and solve forward/backward dynamics problem on it. This package is functioned through:constructing 2d polygon shape rigid bodies and allow motion in 3d space\nconnecting bodies by joints which has 6 degree of freedoms for each\nsolving motions on unconstrained degrees of freedom(passive joints)\nsolving forces on constrained degrees of freedom(active joints)\nplotting/making gif in Julia or making movies in MatlabTo solve a rigid body dynamics problem, this package express the dynamics using 6D spatial vector developed by Roy Featherstone. The governing equations are formulated to fit in half explicit Runge-Kutta method on index-2 differential equations. Constrained forces on joints are represented in Lagrange multiplier terms and solved together with motions of all degrees of freedom.Based on the calculation of dynamical systems, Dyn3d.jl can also be used to simulate fluid-structure interaction problems using package Whirl.jl, with strongly coupled method(finished) and fully coupled method(in package FSI.jl).(Image: )"
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "This package requires Julia 0.6. To install, simply download this Github repository, find the location of this repository expressed in Julia byjulia> Dyn3d_dir = Pkg.dir(\"Dyn3d\")and then setup a symbolic link in shell followingshell$ sudo ln -s actual_address Dyn3d_dirThe plots in this documentation are generated using Plots.jl. You might want to install that too to follow the examples.Or a simple solution asinclude(path*\"Dyn3d.jl\")"
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
    "text": ""
},

{
    "location": "manual/construct_system.html#User-oriented-set-up-1",
    "page": "Construct body-joint system",
    "title": "User-oriented set up",
    "category": "section",
    "text": "ConfigDataType"
},

{
    "location": "manual/construct_system.html#Code-oriented-construction-1",
    "page": "Construct body-joint system",
    "title": "Code-oriented construction",
    "category": "section",
    "text": "Before running any dynamics, we need to construct the body-joint system, connecting bodies by joints in the correct hierarchy using supplied configuration information. In order to do that, three key functions are used:AddBody(id::Int, cb::ConfigBody)\nAddJoint(id::Int, cj::ConfigJoint)\nAssembleSystem(bs, js, sys)which is responsible for add a body to the body system, add a joint to the joint system and assemble the body-joint system and fill in extra hierarchy information. These three functions are also bundled in the correct sequential order inBuildChain(cbs, cjs, csys)In order to solve for the rigid body-joint chain system, we choose body velocity v and joint displacement qJ as main variables in the time marching scheme. We also need to construct a structure that stores and updates the body-joint chain\'s intermediate variables such as all kinds of transformation matrices X, body position in the inertial space x_i and so on. All information of the body-joint system is bundled into a BodyDyn structure."
},

{
    "location": "manual/construct_system.html#Dyn3d.ConfigDataType.ConfigBody-Tuple{Any}",
    "page": "Construct body-joint system",
    "title": "Dyn3d.ConfigDataType.ConfigBody",
    "category": "method",
    "text": "the final plotting direction is:          ^(y)          |_____>(x)      (-z) the coordinate for verts in 2d input is in [z,x]. So y direction only allow zero-width body\n\n\n\n"
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
    "text": "CutOut2d(bd::BodyDyn,bgs::Vector{BodyGrid})\n\nThis function need to be called only once after GenerateBodyGrid for 2d case of flat plates.\n\nIn Dyn3d, bodies are constructed by quadrilateral/triangles(not lines) in x-z plane for both 2d/3d cases. In Whirl, fluid in 2d cases are constructed in x-y plane. Thus to describe plates as lines in x-y space, we cut out the info on the other sides of the plate. Note that verts are formulated in clockwise direction, with the left-bottom corner as origin.\n\n\n\n"
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
