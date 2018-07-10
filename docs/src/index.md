# Dyn3d

*A 2d/3d rigid body dynamics solver*

The main goal of this repository is to construct a rigid body-joint system
and solve forward/backward dynamics problem on it. This package is functioned
through

- constructing 2d/3d rigid body-joints system
- solving motions on unconstrained degrees of freedom(passive joints) and
forces on constrained degrees of freedom(active joints).
- plotting in Julia and making movies in Matlab

To solve a rigid body dynamics problem, this package express the dynamics using
6D spatial vector developed by Roy Featherstone. The governing equations are
formulated to fit in half explicit Runge-Kutta method on index-2 differential
equations. Constrained forces on joints are represented in Lagrange multiplier
terms and solved together with motions of all degrees of freedom.
