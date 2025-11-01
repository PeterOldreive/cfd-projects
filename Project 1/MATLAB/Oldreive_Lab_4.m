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
    ' mm\nNx = ', num2str(Nx), ' and Ny = ', num2str(Ny), '\nNtot = ', num2str(Ntot), '\n'])

%% Define matrix and vectors

A = zeros(Ntot, Ntot); % Matrix of coefficents
T = zeros(1, Ntot); % Vector of unknowns 
B = zeros(Ntot, 1); % Vector of knowns

%% Fill matrix of coefficents and vector of knowns

for j = 1:Ny % For all y points
    for i = 1:Nx % For all x points 
         % Define coefficents based on discretized equations
         % Corner points
         if(i==1 && j==1) % Bottom left corner
            A(j,i) = 1; 
            A(j,i+1) = -1/2;
            A(j+1, i) = -1/2;
        elseif(i == Nx && j == 1) % Bottom right corner
            A(j,i) = 1; 
            A(j,i-1) = -1/4;
            A(j+1,i) = -1/4;
            B(j) = 1100/4; 
        elseif(i == 1 && j == Ny) % top left corner
            A(j,i) = 1; 
            A(j,i+1) = -1/3;
            A(j-1, i) = -1/3;
            B(j) = 400/3; 
         elseif(i == Nx && j == Ny) % top right corner
            A(i,j) = 1; 
            A(i-1,j) = -1/4;
            A(i,j-1) = -1/4;
            B(j) = 1100/3; 
         end 

         % Boundary nodes
        % if(i*deltax <= xdim/2 && j = 1)



    end
end