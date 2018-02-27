% This file plots body chain position in figure or movie
clc;clear all;

% problem definition
system.ndim = 2;
system.nbody = 2;
for i = 1:system.nbody
    system.body(i).nvert = 4;
end

% construct system structure from input .dat file
system.data=load('verts.dat');

% plot using artic_movie_3d
artic_movie_3d(system,2,'savemovie','movie.avi');