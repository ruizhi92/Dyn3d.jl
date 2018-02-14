cmake_minimum_required(VERSION 2.6)
project(DYN3D)
enable_language(Fortran)
# Ruizhi Yang, 2017 Aug


#---------------------------------------------------------------
# Set flags
#---------------------------------------------------------------

# optimized flag
#set (CMAKE_Fortran_FLAGS "-O2 -cpp -w -ffree-form -x f95-cpp-input -ffree-line-length-none")

#debug flag
set (CMAKE_Fortran_FLAGS "-g -ffree-form -O2  -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace")

#-ffpe-trap=zero,overflow,underflow

#set (CMAKE_Fortran_FLAGS "-g -ffree-form -Wall -fbounds-check -fcheck=all")

#---------------------------------------------------------------
# Set source and subdirectory
#---------------------------------------------------------------

file(GLOB_RECURSE DYN3D_SRC
    "dyn3d.f90"
    "config_*.f90"
    "module_*.f90"
)

add_subdirectory(basic_matrix_operations)
add_subdirectory(ode_methods)

#---------------------------------------------------------------
# build executable
#---------------------------------------------------------------

add_executable(dyn3d ${DYN3D_SRC})
SET_TARGET_PROPERTIES(dyn3d PROPERTIES LINKER_LANGUAGE Fortran Fortran_FORMAT FREE)
