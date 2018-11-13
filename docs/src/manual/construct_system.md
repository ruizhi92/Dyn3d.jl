# Construct body-joint system

## User-oriented set up
`ConfigDataType`

## Code-oriented construction

Before running any dynamics, we need to construct the body-joint system, connecting
bodies by joints in the correct hierarchy using supplied configuration information.
In order to do that, three key functions are used:
- [`AddBody(id::Int, cb::ConfigBody)`](@ref AddBody)
- [`AddJoint(id::Int, cj::ConfigJoint)`](@ref AddJoint)
- [`AssembleSystem(bs, js, sys)`](@ref AssembleSystem)

which is responsible for add a body to the body system, add a joint to the joint system and
assemble the body-joint system and fill in extra hierarchy information. These three functions
are also bundled in the correct sequential order in
- [`BuildChain(cbs, cjs, csys)`](@ref BuildChain)

In order to solve for the rigid body-joint chain system, we choose body velocity $v$
and joint displacement $qJ$ as main variables in the time marching scheme. We also
need to construct a structure that stores and updates the body-joint chain's intermediate
variables such as all kinds of transformation matrices $X$, body position in the inertial
space $x_i$ and so on. All information of the body-joint system is bundled into a `BodyDyn`
structure.








## Methods
```@autodocs
Modules = [ConfigDataType, ConstructSystem]
Order   = [:type, :function]
```

## Index
```@index
Pages = ["construct_system.md"]
```
