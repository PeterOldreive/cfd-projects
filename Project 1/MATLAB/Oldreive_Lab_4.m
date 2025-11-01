clear; close all; clc;
%% 
% Peter Oldreive
% B00894035
% MECH 4865 - CFD
% Project 1
% Created: Nov 1, 2025
% Updated: 
% 
% This program calculates the steady state temperature field of a 2D plane 
% subject to known boundary conditions. This is to statify the requirements 
% of MECH 4865 Project 1. 

%% Define Dirichlet Boundary Conditions

T_top = 400; % Top Dirichlet condition (K)
T_bottom = 400; % Bottom Dirichlet condition (K)
T_right = 700; % Right Dirichlet condition (K)
T_left = 0;        % Left Dirichlet condtition (K)
cNtoD = 200; % Transition x-value of bottom BC from Neumann to Dirichlet (mm)

%% Define Plane Dimensions

xdim = 400; % x-dimension of plate (mm)
ydim = 200; % y-deimension of plate (mm)

%% Define discetized grid
deltax = input("Input x-grid spacing (mm): "); % Discretized grid x-dimension
deltay = input('Input y-grid spacing (mm): '); % Discretized grid y-dimension

% Find number of grid points excluding those on boundary condition
Nx = xdim/deltax - 1; % Number of points in x-direction
Ny = ydim/deltay - 1; % Number of points in y-direction

Ntot = Nx*Ny; % Total number of points

% Display discretized points to user
fprintf(['For deltax = ', num2str(deltax), ' mm and deltay = ', num2str(deltay), ...
    ' mm\nNx = ', num2str(Nx), ' and Ny = ', num2str(Ny), '\nNtot = ', num2str(Ntot)])

%% 