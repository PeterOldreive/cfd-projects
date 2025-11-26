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
    x_values(x) = (x) * deltax; 
end

% Generate vector of discretized y-values
y_values = zeros(Ny,1);
for y = 1:Ny
    y_values(y) = (y) * deltay; 
end

% Make a grid of coordinates
[X, Y] = meshgrid(x_values, y_values);


% % Generate a vector of bottom boundary, and 50 mm steps to top
% T_bottom_values = T_plot(1, :); % Make botttom boundary a vector
% T_50mm_above_bottom = T_plot(fix(50/deltay), :); % Closest point to 
%                                                   % 50mm above bottom
%                                                   % boundary
% T_100mm_above_bottom = T_plot(fix(100/deltay), :); % Closest point to 
%                                                   % 100 mm above bottom
%                                                   % boundary
% T_150mm_above_bottom = T_plot(fix(150/deltay), :); % Closest point to 
%                                                   % 150mm above bottom
%                                                   % boundary
% T_top_values = T_plot(Ny, :); % Top boundary 

%% Plotting

% % 3D Plot
% figure;
% surf(T_plot,'EdgeColor','k');
% xlabel('x-index');
% ylabel('y-index');
% zlabel('Temperature');
% title('Temperature field (surface view)');
% colorbar;
% view(45,45);

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

% % Temperature along bottom boundary, plus 50 mm above 
% figure(); 
% plot(x_values, T_bottom_values, 'o', 'MarkerSize', 6, 'color', 'r');
% hold on;
% plot(x_values, T_50mm_above_bottom, 'o', 'MarkerSize', 6, 'color', 'b');
% hold on; 
% plot(x_values, T_100mm_above_bottom, 'o', 'MarkerSize', 6, 'color', 'g');
% hold on;
% plot(x_values, T_150mm_above_bottom, 'o', 'MarkerSize', 6, 'color', 'k');
% hold on;
% plot(x_values, T_top_values, 'o', 'MarkerSize', 6, 'color', 'm');
% hold on;
% xlabel('x (mm)');
% ylabel('Temperature (°C)');
% title('Temperature Distribution along Bottom Boundary Nodes and 50 mm Above');
% legend('Bottom Boundary', 'y = 50 mm', 'y = 100 mm', 'y = 150 mm', ... 
%     'Top Boundary', 'Location','southoutside'); 
% grid on;
