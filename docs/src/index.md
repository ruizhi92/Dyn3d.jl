# Dyn3d.jl

*A 2d/3d rigid body dynamics solver*

The main goal of this repository is to construct a rigid body-joint system
and solve forward/backward dynamics problem on it. This package is functioned
through:

- constructing 2d polygon shape rigid bodies and allow motion in 2d/3d space
- connecting bodies by joints which has 6 degree of freedoms for each
- solving motions on unconstrained(passive) degrees of freedom of joints
- solving forces on constrained degrees of freedom, allowing active motion
- plotting/making gif in *Julia* or making movies in *Matlab*

To solve a rigid body dynamics problem, this package express the dynamics using
6D spatial vector developed by Roy Featherstone[^1]. The governing equations are
formulated to fit in half explicit Runge-Kutta method on index-2 differential
equations[^2]. Constrained forces on joints are represented in Lagrange multiplier
terms and solved together with motions of all degrees of freedom.

Based on the calculation of dynamical systems, `Dyn3d.jl` is also used to simulate
fluid-structure interaction(FSI) problems together with package `Whirl.jl` for
strongly coupled method. Notebook example is provided in notebook folder. Fully
coupled method taking advantage of both `Dyn3d.jl` and `Whirl.jl` is implemented
in package `FSI.jl`.

![](https://github.com/ruizhi92/Dyn3d.jl/raw/master/example_gif.gif)

## Installation

This package currently requires *Julia* 0.6.
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

If you have trouble in setting up the symbolic like to directory of `Dyn3d.jl`,
a simple alternative solution is:
```
include(path*"Dyn3d.jl")
```

## References

[^1]: Featherstone, Roy. Rigid body dynamics algorithms. Springer, 2014.
[^2]: Brasey, Valérie, and Ernst Hairer. "Half-explicit Runge–Kutta methods for differential-algebraic systems of index 2." SIAM Journal on Numerical Analysis 30, no. 2 (1993): 538-552.
