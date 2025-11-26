clear all; close all; clc; tic;
%% 
% Peter Oldreive
% B00894035
% MECH 4865 - CFD
% Project 2
% Created: Nov 24, 2025
% Updated: Nov 26, 2025
% 
% This program calculates the steady state temperature field of a 2D plane 
% subject to known boundary conditions. This is to statify the requirements 
% of MECH 4865 Project 2. The program does not use a direct solver to
% calcualte the temperature distribution, the temperatures are calculated
% using the Gauss_Seidel method by referencing an external script. 

%% Define Plate dimensions 
xdim = 800; % x-dimension of plate (mm)
ydim = 400; % y-deimension of plate (mm)           

%% Define Boundary Conditions

% Dirichlet Boundaries 
T_top = 0; % Top Dirichlet condition (K)
T_bottom = 0; % Bottom Dirichlet condition (K)
T_right = 0; % Right Dirichlet condition (K)
T_left = 0;        % Left Dirichlet condtition (K)

% Neumann conditions
% West Wall
% South wall from 0 <= x <= 400 mm
% East Wall
bc_trasnition = 400; % Transition point of bottom boundary condtion from 
                     % Neumann condtion to Newton/Robin

% Newton/Robin Boundaries
h_south = 1200; % South boundary convection coefficient (W/m^2)
T_inf_south = 300; % South boundary fluid temperature (K) 

h_north = 1500; % North boundary convection coefficent (W/m^2)
T_inf_north = 1800; % North boundary fluid temperature (K) 

%% Define Finite Volume Cell Size 
while(1)
    deltax = input("Cell x-dimension (mm): "); % Cell x-dimension (mm)
    deltay = input("Cell y-dimension (mm): "); % Cell y-dimension (mm)
    
    % Check for valid grid spacing
    if(mod(xdim, deltax) == 0 && mod(ydim, deltay) == 0)
        % Ensure the dimensions are divisable by the given cell dimensions
        break; % Exit the loop if valid grid spacing is provided
    else
        disp('Invalid grid spacing');
    end
end
% Define number of cells in the grid
Nx = xdim/deltax; % Number of cells in x-dimension 
Ny = ydim/deltay; % Number of cells in y-dimension 

Ntot = Nx*Ny; % Total Number of cells 

% Display Grid information to user 
fprintf(['For deltax = ', num2str(deltax), ' mm and deltay = ', num2str(deltay), ...
    ' mm\nNx = ', num2str(Nx), ' and Ny = ', num2str(Ny), '\nNtot = ', num2str(Ntot), '\n'])

%% Define matrix and vectors

A = zeros(Ntot, Ntot); % Matrix of coefficents
T = zeros(Ntot, 1); % Vector of unknowns 
B = zeros(Ntot, 1); % Vector of knowns

%% Fill matrix of coefficents and vector of knowns

% NEED TO FILL OUT %

for i = 1:Ntot % For all points 
j = fix((i-1)/Nx) + 1; % Calculate the row number
  % Corner points 
    if(i == 1) % Bottom left corner
       A(i, i) = 2; % Point
       A(i, i+1) = -1; % East point
       A(i, i+Nx) = -1; % North point
    elseif(i == Nx) % Bottom right corner
       A(i, i) = 4; % Point
       A(i, i-1) = -1; % West point
       A(i, i+Nx) = -1; % North point
       B(i) = (T_right + T_bottom); % Known
    elseif(i == Ntot - Nx + 1) % Top Left
        A(i, i) = 3; % Point
        A(i, i+1) = -1; % East Point
        A(i, i-Nx) = -1; % South Point
        B(i) = (T_top); % Known
    elseif(i == Ntot) % Top right
        A(i, i) = 4; % Point
        A(i, i-1) = -1; % West point
        A(i, i-Nx) = -1; % South Point
        B(i) = (T_top + T_right); % Known
    % Boundary Nodes
    elseif(j == 1 && i*deltax < xdim/2 && i > 1) % Left Lower boundary
        A(i, i) = 3; % Point
        A(i, i+1) = -1; % East point
        A(i, i-1) = -1; % West point
        A(i, i+Nx) = -1; % North point
    elseif(j == 1 && i*deltax >= xdim/2 && i < Nx) %Right lower boundary
        % Right lower boundary contains node at transition point
        %from Neumann condition to Dirichlet condition. 
        A(i, i) = 4; % Point
        A(i, i+1) = -1; % East point
        A(i, i-1) = -1; % West point
        A(i, i+Nx) = -1; % North point
        B(i) = T_bottom; % Known
    elseif(j > 1 && mod(i, Nx) == 1 && j<Ny) % Left boundary
        A(i, i) = 3; % Point
        A(i, i+1) = -1; % East point
        A(i, i+Nx) = -1; % North point
        A(i, i-Nx) = -1; % South point
    elseif(j > 1 && mod(i, Nx) == 0 && j < Ny) % Right boundary
        A(i, i) = 4; % Point
        A(i, i-1) = -1; % West point
        A(i, i+Nx) = -1; % North point
        A(i, i-Nx) = -1; % South point
        B(i) = T_right; % Known
    elseif(j == Ny && i > Nx*(j-1) + 1 && i < Nx*j) % Top 
        A(i, i) = 4; % Point
        A(i, i-1) = -1; % West point
        A(i, i+1) = -1; % East point
        A(i, i-Nx) = -1; % South point
        B(i) = T_top; % Known
    else % Interior nodes
        A(i, i) = 4; % Point
        A(i, i-1) = -1; % West point
        A(i, i+1) = -1; % East point
        A(i, i+Nx) = -1; % North point
        A(i, i-Nx) = -1; % South point
    end
end

%% Call Solver Script 

writematrix(A, "tools/coeff_matrix.csv"); % Write Matrix of Coefficients (A) 
                                          % to a .csv file
writematrix(B, 'tools/knowns.csv'); % Write vector of knowns (B) to a .csv file

% Run solver script saved one directory below within the tools folder.
% changing directory may be required to run script on a different PC
run('tools/gauss_seidel.m');

% Read in solution vector gererated by gauss seidel mechanism 
T = readmatrix('tools/solution_vector.csv'); 

%Record compute time
endtime = toc; 
fprintf("Compute time: %d s", endtime) % Print compute time 

%% Post Process 
