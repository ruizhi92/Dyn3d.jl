This is the 3d dynamic chain solver using 6d spatial vector.

Code is written in Fortran.

branch **master** computes dynamics using articulated body method.

branch **HERK_local_coord** computes dynamics formulating into a half-explicit Runge Kutta method solver.

branch **HERK_body** tried to use a body-orientated approach but didn't finish because of some fuzzy concept.

branch **HERK_inertial_coord** tried to solve the system in an absolute inertial coodinate, but didn't work.