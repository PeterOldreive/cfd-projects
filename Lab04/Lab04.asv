clc; clear all; close all;


Nx = 15;
Ny = 15; 

N = NX*Ny;

A = zeros(N, N);
B = zerps(N, 1);

for j = 1:Ny
    for i = 1:Nx 
        k = index(i.j); 
        A(k, k) = -4; 

        % East Neighbor
        if i < Nx
            A(k, )