# Dyn3d.jl

*A 2d/3d rigid body dynamics solver*

The main goal of this repository is to construct a rigid body-joint system
and solve forward/backward dynamics problem on it. This package is functioned
through:

- constructing 2d polygon shape rigid bodies and allow motion in 3d space
- connecting bodies by joints which has 6 degree of freedoms for each
- solving motions on unconstrained degrees of freedom(passive joints)
- solving forces on constrained degrees of freedom(active joints)
- plotting/making gif in *Julia* or making movies in *Matlab*

To solve a rigid body dynamics problem, this package express the dynamics using
6D spatial vector developed by Roy Featherstone. The governing equations are
formulated to fit in half explicit Runge-Kutta method on index-2 differential
equations. Constrained forces on joints are represented in Lagrange multiplier
terms and solved together with motions of all degrees of freedom.

Based on the calculation of dynamical systems, `Dyn3d.jl` can also be used to simulate
fluid-structure interaction problems using package `Whirl.jl`, with strongly coupled
method(finished) and fully coupled method(in package `FSI.jl`).

![](https://github.com/ruizhi92/Dyn3d.jl/raw/master/example_gif.gif)

## Installation

This package requires *Julia* 0.6.
To install, simply download this Github repository, find the location of this repository
expressed in *Julia* by
```
julia> Dyn3d_dir = Pkg.dir("Dyn3d")
```
and then setup a symbolic link in shell following
```
shell$ sudo ln -s actual_address Dyn3d_dir
```

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).
You might want to install that too to follow the examples.

Or a simple solution as
```
include(path*"Dyn3d.jl")
```
