# Dyn3d

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://ruizhi92.github.io/Dyn3d.jl/latest)
[![Build Status](https://travis-ci.org/ruizhi92/Dyn3d.jl.png?branch=master)](https://travis-ci.org/ruizhi92/Dyn3d.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/0ykpxm3e8rftro6m/branch/master?svg=true)](https://ci.appveyor.com/project/ruizhi92/dyn3d-jl/branch/master)
[![codecov](https://codecov.io/gh/ruizhi92/Dyn3d.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ruizhi92/Dyn3d.jl)

## About the package

This is the 2d/3d rigid body dynamics solver using 6d spatial vector. Examples notebooks
are given in notebook folder. User just need to change the configuration files
for different cases, nothing else needed to be changed.

Code is written in Julia and Fortran on different branch.

- branch **master** for `Julia 1.1`
- branch **v0.6** for `Julia 0.6`
- branch **v0.7** for `Julia 0.7`
- branch **Fortran/artic3d** computes dynamics using articulated body method.
- branch **Fortran/HERK** computes dynamics formulating into a half-explicit Runge Kutta method solver in Fortran.

**Dyn3d.jl** is registered in the general Julia registry. To install, type
e.g.,
```julia
] add Dyn3d
```

Then, in any version, type
```julia
julia> using Dyn3d
```
See the example Jupyter notebooks in the examples folder.

![](https://github.com/ruizhi92/Dyn3d.jl/raw/master/example_gif.gif)
