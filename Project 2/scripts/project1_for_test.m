clear; close all; clc; tic;
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
T = zeros(Ntot, 1); % Vector of unknowns 
B = zeros(Ntot, 1); % Vector of knowns

%% Fill matrix of coefficents and vector of knowns

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

%%  Direct solve for T values
T = A \ B; 
  
csvwrite("tools/coeff_matrix.csv", A);
csvwrite('tools/knowns.csv', B);

run('tools/gauss_seidel.m');

TGS = readmatrix('tools/solution_vector.csv');

Test = abs(TGS - T); 
max(Test)

%Record compute time
endtime = toc; 
fprintf("Compute time: %d s \n", round(endtime)) % Print compute time 

%% Post Process

% Turn vector T into an NyxNx matrix
T_plot = zeros(Ny, Nx);
index = 1; 
for y = 1:Ny 
    for x = 1:Nx
     T_plot(y,x) = T(index); 
     index = index + 1;
    end
end
writematrix(T_plot, 'Temperature_Field.csv');

% Generate vector of discretized x-values
x_values = zeros(Nx,1);
for x = 1:Nx
    x_values(x) = (x) * deltax; 
end

% Generate vector of discretized y-values
y_values = zeros(Ny,1);
for y = 1:Ny
    y_values(y) = (y) * deltay; 
end

% Make a grid of coordinates
[X, Y] = meshgrid(x_values, y_values);

% Generate a vector of bottom boundary, and 50 mm steps to top
T_bottom_values = T_plot(1, :); % Make botttom boundary a vector
T_50mm_above_bottom = T_plot(fix(50/deltay), :); % Closest point to 
                                                  % 50mm above bottom
                                                  % boundary
T_100mm_above_bottom = T_plot(fix(100/deltay), :); % Closest point to 
                                                  % 100 mm above bottom
                                                  % boundary
T_150mm_above_bottom = T_plot(fix(150/deltay), :); % Closest point to 
                                                  % 150mm above bottom
                                                  % boundary
T_top_values = T_plot(Ny, :); % Top boundary 


%% Plotting

% 3D Plot
figure;
surf(T_plot,'EdgeColor','k');
xlabel('x-index');
ylabel('y-index');
zlabel('Temperature');
title('Temperature field (surface view)');
colorbar;
view(45,45);

% Show discretized 2D grid of Nx x Ny points
figure;
plot(X, Y, 'o', 'LineWidth', 2, 'color', 'r', 'MarkerSize',3)
set(gca,'YDir','normal'); % make y increase upward visually
axis equal tight;
xlim([0 xdim])
ylim([0 ydim])
title('Discretized Grid of Nx x Ny Points Spaced by \Deltax and \Deltay');
xlabel('x-Coordinate (mm)');
ylabel('y-Coordinate (mm)'); 


% 2D Temperature Field
figure;
imagesc(T_plot);
hold on
set(gca,'YDir','normal'); % make y increase upward visually
axis equal tight;
colorbar;
title('Temperature contour map');
xlabel('x-Coordinate (mm)');
ylabel('y-Coordinate (mm)'); 
% plot(X/10, Y/10, 'o', 'LineWidth', 2, 'color', 'r', 'MarkerSize',3)

% Temperature along bottom boundary, plus 50 mm above 
figure(); 
plot(x_values, T_bottom_values, 'o', 'MarkerSize', 6, 'color', 'r');
hold on;
plot(x_values, T_50mm_above_bottom, 'o', 'MarkerSize', 6, 'color', 'b');
hold on; 
plot(x_values, T_100mm_above_bottom, 'o', 'MarkerSize', 6, 'color', 'g');
hold on;
plot(x_values, T_150mm_above_bottom, 'o', 'MarkerSize', 6, 'color', 'k');
hold on;
plot(x_values, T_top_values, 'o', 'MarkerSize', 6, 'color', 'm');
hold on;
xlabel('x (mm)');
ylabel('Temperature (°C)');
title('Temperature Distribution along Bottom Boundary Nodes and 50 mm Above');
legend('Bottom Boundary', 'y = 50 mm', 'y = 100 mm', 'y = 150 mm', ... 
    'Top Boundary', 'Location','southoutside'); 
grid on;

