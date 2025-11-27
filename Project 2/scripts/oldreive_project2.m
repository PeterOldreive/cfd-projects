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
xdim = 0.800; % x-dimension of plate (m)
ydim = 0.400; % y-deimension of plate (m)           

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

% Thermal coductivity of the plate 
k = 80; % Thermal conductivity of the plate (W/mK)

%% Define Finite Volume Cell Size 
while(1)
    deltax = input("Cell x-dimension (mm): ")/1000; % Cell x-dimension (m)
    deltay = input("Cell y-dimension (mm): ")/1000; % Cell y-dimension (m)

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
yx = deltay/deltax; % Shorthand coefficient for use later
xy = deltax/deltay; % Shorthand coefficient for use later 
% Define coefficents from ficticious node analysis 

% North boundary
apfn = 1 + (deltay*h_north)/(2*k); % coefficent from equation derivation
asfn = -1 + (deltay*h_north)/(2*k); % coefficent from equation derivation
bpfn = (h_north*deltay*T_inf_north)/k; % coefficent from equation derivation

% South boundary 
apfs = 1 + (deltay*h_south)/(2*k); % coefficent from equation derivation
anfs = - 1 + (deltay*h_south)/(2*k); % coefficent from equation derivation
bpfs = (h_south*deltay*T_inf_south)/k; % coefficent from equation derivation

for i = 1:Ntot % For all points 
j = fix((i-1)/Nx) + 1; % Calculate the row number
  % Corner points 
    if(i == 1) % Bottom left corner
       A(i, i) = -(xy + yx); % Cell
       A(i, i+1) = yx; % East Cell
       A(i, i+Nx) = xy; % North cell 
    elseif(i == Nx) % Bottom right corner
       A(i, i) = -((2+anfs/apfs)*xy + yx); % Cell
       A(i, i-1) = yx; % West Cell
       A(i, i+Nx) = xy; % North Cell
       B(i) = -xy*bpfs/apfs; % Known
    elseif(i == Ntot - Nx + 1) % Top Left
        A(i, i) = -((2+asfn/apfs)*xy + yx); % Cell
        A(i, i+1) = yx; % East Cell
        A(i, i-Nx) = xy; % South Cell
        B(i) = -xy*bpfn/apfn; % Known
    elseif(i == Ntot) % Top right
        A(i, i) = -((2+asfn/apfn)*xy + yx); % Point
        A(i, i-1) = yx; % West point
        A(i, i-Nx) = xy; % South Point
        B(i) = -xy*bpfn/apfn; % Known
    % Boundary Nodes
    elseif(j == 1 && i*deltax < xdim/2 && i > 1) % Left Lower boundary
        A(i, i) = -(2*yx+xy); % Point
        A(i, i+1) = yx; % East point
        A(i, i-1) = xy; % West point
        A(i, i+Nx) = xy; % North point
    elseif(j == 1 && i*deltax >= xdim/2 && i < Nx) %Right lower boundary
        % Right lower boundary contains node at transition point
        %from Neumann condition to Newton/Robin Condition. 
        A(i, i) = -((2+anfs/apfs)*xy + 2*yx); % Point
        A(i, i+1) = yx; % East point
        A(i, i-1) = yx; % West point
        A(i, i+Nx) = xy; % North point
        B(i) = -xy*bpfs/apfs; % Known
    elseif(j > 1 && mod(i, Nx) == 1 && j<Ny) % Left boundary
        A(i, i) = -(2*xy + yx); % Point
        A(i, i+1) = yx; % East point
        A(i, i+Nx) = xy; % North point
        A(i, i-Nx) = xy; % South point
    elseif(j > 1 && mod(i, Nx) == 0 && j < Ny) % Right boundary
        A(i, i) = -(2*xy + yx); % Point
        A(i, i-1) = yx; % West point
        A(i, i+Nx) = xy; % North point
        A(i, i-Nx) = xy; % South point
    elseif(j == Ny && i > Nx*(j-1) + 1 && i < Nx*j) % Top 
        A(i, i) = -((2+asfn/apfn)*xy + 2*yx); % Point
        A(i, i-1) = yx; % West point
        A(i, i+1) = yx; % East point
        A(i, i-Nx) = xy; % South point
        B(i) = -xy*bpfn/apfn; % Known
    else % Interior nodes
        A(i, i) = -2*(xy+yx); % Point
        A(i, i-1) = yx; % West point
        A(i, i+1) = yx; % East point
        A(i, i+Nx) = xy; % North point
        A(i, i-Nx) = xy; % South point
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
fprintf("Compute time: %d s\n", round(endtime)) % Print compute time





%% Test solution
T2 = A \ B; 
Test = abs(T2 - T); 
fprintf('%d \n', max(Test))

%% Post Process 

% Turn temperature vector into a NyxNx matrix for plotting

T_plot = zeros(Ny, Nx); % Define matrix
index = 1; % Index prepresents point number  
for y = 1:Ny % For all rows
    for x = 1:Nx % For all columns
     T_plot(y,x) = T(index);  
     index = index + 1;
    end
end

% Write matrix of temperature field for later comparison
writematrix(T_plot, 'Temperature_Field.csv');

% Generate vector of discretized x-values
x_values = zeros(Nx,1);
for x = 1:Nx 
    x_values(x) = (x-1) * deltax; 
end

% Generate vector of discretized y-values
y_values = zeros(Ny,1);
for y = 1:Ny
    y_values(y) = (y-1) * deltay; 
end

% Make a grid of coordinates
[X, Y] = meshgrid(x_values, y_values);


% Generate a vector of bottom boundary, and 100 mm steps to top
T_bottom_values = T_plot(1, :); % Make botttom boundary a vector
T_100mm_above_bottom = T_plot(fix(0.100/deltay), :); % Closest point to 
                                                  % 50mm above bottom
                                                  % boundary
T_200mm_above_bottom = T_plot(fix(0.200/deltay), :); % Closest point to 
                                                  % 100 mm above bottom
                                                  % boundary
T_300mm_above_bottom = T_plot(fix(0.300/deltay), :); % Closest point to 
                                                  % 150mm above bottom
                                                  % boundary
T_top_values = T_plot(Ny, :); % Top boundary 

% Generate cutlines paralell to the y-axis in 100 mm increments
T_left_values = T_plot(:,1); % 1st column is left boundary
T_x_200 = T_plot(:, fix(0.200/deltax)); % Temperatures at 100 mm right of the
                                      % left boundary
T_x_400 = T_plot(:, fix(0.400/deltay)); % temperatures at 200 mm right of the 
                                      % left boundary
T_x_600 = T_plot(:, fix(0.600/deltay)); % Temperatures at 300 mm right of the 
                                      % left boundary
T_right_values = T_plot(:, Nx); % Right boundary temperatues

% Find diagonal vector
rowcount = 1; % Start in 1st row
T_diag = zeros(Nx/2, 1); % Define diagonal vector
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
xlabel('x-Coordinate (m)');
ylabel('y-Coordinate (m)'); 

% 2D Colourmap
figure('Name', '2D Colour Map','NumberTitle','off');
hold on
contourf(T_plot,500:20:1800,'LineColor', "none"); % Add Contour Lines every 10 K
set(gca,'YDir','normal');
axis equal tight;
a=colorbar;
a.Label.String = 'Temperature, T (K)';
a.Limits = [500 1800];
title('2D Colourmap with Contour Lines');
xlabel('x-Coordinate (m)');
ylabel('y-Coordiante (m)');
xticks(0:0.10/deltax:Nx)
xticklabels(0:0.10:xdim)
yticks(0:0.1/deltay:Ny)
yticklabels(0:0.10:ydim)


% Temperature along horizontal cutlines 
figure(); 
plot(x_values, T_bottom_values, 'o', 'MarkerSize', 6, 'color', 'r', ...
    'LineStyle','--' );
hold on;
plot(x_values, T_100mm_above_bottom, 'o', 'MarkerSize', 6, 'color', 'b', ...
    'LineStyle','--' );
hold on; 
plot(x_values, T_200mm_above_bottom, 'o', 'MarkerSize', 6, 'color', 'g', ...
    'LineStyle','--' );
hold on;
plot(x_values, T_300mm_above_bottom, 'o', 'MarkerSize', 6, 'color', 'k', ...
    'LineStyle','--' );
hold on;
plot(x_values, T_top_values, 'o', 'MarkerSize', 6, 'color', 'm', ...
    'LineStyle','--');
hold on;
xlabel('x (m)');
ylabel('Temperature (K)');
title('Temperature Distribution along Horizontal Cutlines');
legend('Bottom Boundary', 'y = 100 mm', 'y = 200 mm', 'y = 300 mm', ... 
    'Top Boundary', 'Location','southwest'); 
grid on;

% Temperature cutlines paralell to y-axis
figure(); 
plot( y_values, T_left_values, 'o', 'MarkerSize', 6, 'color', 'r', ...
    'LineStyle','--');
hold on;
plot( y_values, T_x_200, 'o', 'MarkerSize', 6, 'color', 'b', ...
    'LineStyle','--');
hold on; 
plot(y_values, T_x_400, 'o', 'MarkerSize', 6, 'color', 'g', ... 
    'LineStyle','--');
hold on;
plot(y_values, T_x_600, 'o', 'MarkerSize', 6, 'color', 'k', ...
    'LineStyle', '--');
hold on;
plot(y_values , T_right_values, 'o', 'MarkerSize', 6, 'color', 'm', ...
    'LineStyle', '--');
hold on;
xlabel('y-Coordinate (m)');
ylabel('Temperature (K)');
title('Temperature Distribution along Veritical Cutlines');
legend('Left Boundary', 'x = 200 mm', 'x = 400 mm', 'x = 600 mm', ... 
    'Right Boundary', 'Location','southeast'); 
grid on;

% Diagonal Temperature 
figure(); 
plot(deltax:2*deltax:xdim, T_diag, 'o', 'MarkerSize', 5, 'color', 'r', ...
    'LineStyle','--');
hold on;
ylabel('Temperature (K)');
xlabel('x-Coordinate (mm)');
title('Temperature Along Diagonal Cutline (Bottom Left to Top Right)');
grid on;