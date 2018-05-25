This is the 3d rigid body dynamic chain solver using 6d spatial vector.

Code is written in Julia and Fortran on different branch.

branch **master** computes dynamics formulating into a half-explicit Runge Kutta method solver.

branch **Fortran/artic3d** computes dynamics using articulated body method.

branch **Fortran/HERK** rewrites **master** in Fortran.
