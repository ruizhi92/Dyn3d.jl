# Dyn3d

[![Build Status](https://travis-ci.org/ruizhi92/Dyn3d.jl.png?branch=master)](https://travis-ci.org/ruizhi92/Dyn3d.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://ruizhi92.github.io/Dyn3d.jl/latest)
[![codecov](https://codecov.io/gh/ruizhi92/Dyn3d.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ruizhi92/Dyn3d.jl)

This is the 3d rigid body dynamic chain solver using 6d spatial vector.

Code is written in Julia and Fortran on different branch.

- branch **master** computes dynamics formulating into a half-explicit Runge Kutta method solver.
- branch **Fortran/artic3d** computes dynamics using articulated body method.
- branch **Fortran/HERK** rewrites **master** in Fortran.

This package's local dir will be set by user. Find Julia repo address by
julia> Pkg.dir("Dyn3d")
Then you can make a symlinking by
shell$ sudo ln -s actual_address Julia_repo_address
