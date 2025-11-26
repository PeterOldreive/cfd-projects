clear; close all; clc;
%% 
% Peter Oldreive
% B00894035
% MECH 4865 - CFD
% Project 1
% Created: Nov 1, 2025
% Updated: Nov 5, 2025
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
% writematrix(T_plot, 'Temperature_Field.csv');

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

% Generate cutlines paralell to the y-axis in 100 mm increments
T_left_values = T_plot(:,1); % 1st column is left boundary
T_x_100 = T_plot(:, fix(100/deltax)); % Temperatures at 100 mm right of the
                                      % left boundary
T_x_200 = T_plot(:, fix(200/deltay)); % temperatures at 200 mm right of the 
                                      % left boundary
T_x_300 = T_plot(:, fix(300/deltay)); % Temperatures at 300 mm right of the 
                                      % left boundary
T_right_values = T_plot(:, Nx); % Right boundary temperatues

% Find diagonal vector
rowcount = 1; % Start in 1st row
for columncount = 2:Nx % count up to Nx
   
    if (mod(columncount, 2) == 0 && rowcount < Ny+2) % count up rows until max
        T_diag(rowcount, 1) = T_plot(rowcount, columncount); % index diagonals out of matrix
        rowcount = rowcount + 1; % Count up row     
    end

  
end 



%% Plotting

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

% 2D Colourmap
figure('Name', '2D Colour Map','NumberTitle','off');
hold on
contourf(T_plot,380:10:650,'LineColor', "none"); % Add Contour Lines every 10 K
set(gca,'YDir','normal');
axis equal tight;
a=colorbar;
a.Label.String = 'Temperature, T (K)';
a.Limits = [400 700];
title('2D Colourmap with Contour Lines');
xlabel('x-Coordinate (mm)');
ylabel('y-Coordiante (mm)');
xticks(0:50/deltax:xdim/deltax)
xticklabels(0:50:xdim)
yticks(0:50/deltay:ydim/deltay)
yticklabels(0:50:ydim)
xlim([0 xdim/deltax]);
ylim([0 ydim/deltay]); 


% Temperature cutlines paralell to x-axis
figure(); 
plot(x_values, T_bottom_values, '.', 'MarkerSize', 8, 'color', 'r', ...
    'LineStyle','--');
hold on;
plot(x_values, T_50mm_above_bottom, '+', 'MarkerSize', 6, 'color', 'b', ...
    'LineStyle','--' );
hold on; 
plot(x_values, T_100mm_above_bottom, '*', 'MarkerSize', 6, 'color', 'g', ...
    'LineStyle', '--');
hold on;
plot(x_values, T_150mm_above_bottom, 'x', 'MarkerSize', 6, 'color', 'k', ...
    'LineStyle','--');
hold on;
plot(x_values, T_top_values, '^', 'MarkerSize', 6, 'color', 'm', ...
    'LineStyle','--');
hold on;
ylim([350 700])
xlabel('x-Coordinate (mm)');
ylabel('Temperature (K)');
title('Temperature along Horizontal Cutlines');
legend('Bottom Boundary', 'y = 50 mm', 'y = 100 mm', 'y = 150 mm', ... 
    'Top Boundary', 'Location','northwest'); 
grid on;


% Temperature along bottom boundary near condition change
figure(); 
plot(x_values, T_bottom_values, '.', 'MarkerSize', 8, 'color', 'r', ...
    'LineStyle','--');
xlim([170 230])
ylim([390 420])
xlabel('x-Coordinate (mm)');
ylabel('Temperature (K)');
title('Temperature along Bottom Boundary Near Condition Change');
grid on;

% Diagonal Temperature 
figure(); 
plot(deltax*2:2*deltax:xdim-deltax, T_diag, 'o', 'MarkerSize', 5, 'color', 'r', ...
    'LineStyle','--');
hold on;
ylabel('Temperature (K)');
xlabel('x-Coordinate (mm)');
title('Temperature Along Diagonal Cutline (Bottom Left to Top Right)');
grid on;

% Temperature cutlines paralell to y-axis
figure(); 
plot( y_values, T_left_values, '.', 'MarkerSize', 8, 'color', 'r', ...
    'LineStyle','--');
hold on;
plot( y_values, T_x_100, '+', 'MarkerSize', 6, 'color', 'b', ...
    'LineStyle','--');
hold on; 
plot(y_values, T_x_200, '*', 'MarkerSize', 6, 'color', 'g', ... 
    'LineStyle','--');
hold on;
plot(y_values, T_x_300, 'x', 'MarkerSize', 6, 'color', 'k', ...
    'LineStyle', '--');
hold on;
plot(y_values , T_right_values, '^', 'MarkerSize', 6, 'color', 'm', ...
    'LineStyle', '--');
hold on;
ylim([350 700])
xlabel('y-Coordinate (mm)');
ylabel('Temperature (K)');
title('Temperature Distribution along Veritical Cutlines');
legend('Left Boundary', 'x = 100 mm', 'x = 200 mm', 'x = 300 mm', ... 
    'Right Boundary', 'Location','southoutside'); 
grid on;


 
