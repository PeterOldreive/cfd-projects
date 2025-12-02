clear all; close all; clc; 
%% 
% Peter Oldreive
% B00894035
% MECH 4865 - CFD
% Project 3
% Created: Dec 1, 2025
% Updated: Dec 2, 2025
% 
% This program calculates the transient response of the oscillating 
% disturbances of a wave function. This program used a finite difference
% Forward Time Backward Space (FTBS) scheme to discretize the governing
% equation and solve with the given initial and boundary conditions 


%% Define Wave function boundary conditions 
function U_BC = U_BC_func(t, x, C) 
    % Check if time is zero, an x is within given boundary
    if(t == 0 && 0 <= x && x <= 1)
        U_BC = sin(2*pi*x); % Boundary at t = 0
    elseif(t >= 0 && x ==0) % Check that x is zero, and time is greater than zero
        U_BC = sin(2*pi*C*t); % Boundary at x = 0 
    end 
end 

function U = calc_U(Uin, Ui_1n, c, deltax, deltat)
    % Function to calculate U, from Uin and Ui-1n, usig FTBD
    % This equation was derived from the wave equation, using Taylor Series
    % expansions to approximate the partial derivative terms with first
    % order accuracy. 
    U = (1-(c*deltat)/deltax)*Uin + ((c*deltat)/deltax)*Ui_1n;
end 

%% Define Finite Volume Cell Size 
%Take user input 
% deltax = input("Grid x-dimension (m): "); % Discretized delta x
% deltat = input("Time step (s): "); % Time step (s)
deltax = 0.01;
deltat = 0.01; 


% Cmin = input('Wave Speed min, Cmin: '); % Minimim Wave speed
% Cmax = input('Wave Speed max, Cmax: '); % Maximum Wave speed
% Cstep = input('Wave Speed Step, Cstep: '); % Wave speed step
Cmin = 0.25; % Minimum Wave Speed  
Cmax = 1; % Maximum Wave Speed
Cstep = 0.25; % Wave speed step


% Intitialize maximum domain values 
xdim = 1; % Maximium x-dimension (m) 
tmax = 4; % Maximum simulation time (s)


% Define Number of discrete x-values
Nx = xdim/deltax; 

% Define number of time steps
Nt = tmax/deltat;

% Display domain information to user 
fprintf(['For deltax = ', num2str(deltax), ' m and deltat = ', num2str(deltat), ...
    ' m\nNx = ', num2str(Nx), ' and Nt = ', num2str(Nt), '\n'])


% Define solution matrix and vectors of independant variables 
U = zeros(Nx,Nt, Cmax/Cstep); % Matrix of Unknowns
t = zeros(Nt, 1); % Vector of timesteps
x = zeros(Nx, 1); % Vector of dscrete x-coordinatres 

% Define boundary conditions 
for C = Cmin:Cstep:Cmax
    for n = 1:Nt % For all time values
        t(n) = (n-1)*deltat; % Calculate time at point n
        U(1,n, C/Cstep) = U_BC_func(t(n), 0, C); % Calculate boundary condition for x = 0, tn
    end 
    for n = 1:Nx % for all discrete x-values 
        x(n) = (n-1)*deltax; % Calculate position at point i 
        U(n,1, C/Cstep) = U_BC_func(0, x(n), C); % Calculate boundary condition for t = 0, xi
    end 

% Fill U matrix 

    for n = 1:Nt -1  % For all time steps
    
        for j = 2:Nx % For all x-coordinates
      
            % Calculate U value using FTBS scheme
            U(j,n + 1, C/Cstep) = calc_U(U(j,n, C/Cstep), U(j-1,n,C/Cstep), C, deltax, deltat);
        end 
    end 
end 

%% Plotting
%% Surface Plot
for C = Cmin:Cstep:Cmax
    % Surface Plot
    figure('Name',['Surface Plot C ' num2str(C)],'NumberTitle','off')
    surf(U(:, :, C/Cstep),'EdgeColor','none')
    %set(gcf, 'Position', [350 250 800 400])
    xlabel('Time');
    ylabel('x-Direction');
    zlabel('Amplitude');
    title('Wave Propogation (Surface View)');
    a=colorbar;
    a.Label.String = 'Amplitude';
    yticks(0:0.1/deltax:xdim/deltax)
    yticklabels(x(1):0.1:xdim);
    xticks(t(1)/deltat:1/deltat:tmax/deltat)
    xticklabels(t(1):1:tmax);
    view(-45,45);
end 

%% Wave Ampltidue WTR time
figure('Name','Time Plots','NumberTitle','off');
xt = tiledlayout(2,2); % 2 rows, 2 columns
xt.TileSpacing = 'compact';
xt.Padding = 'compact';
title(xt, 'Oscillating Displacement WRT Time for Various C and x')
nexttile;
% C = 0.25
plot( t, U(1,:, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
plot( t, U(1,:, 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
plot( t, U(1, :, 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
plot( t, U(1, :, 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold on;
title('x = 0 m');
grid on;
xlabel('Time (s)');
ylabel('Displacement (m)');
ylim([-1 1])

nexttile;
% C = 0.25
p1 = plot( t, U(0.25/deltax,:, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
p2 = plot( t, U(0.25/deltax ,: , 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
p3 = plot( t, U(0.25/deltax, : , 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
p4 = plot( t, U(0.25/deltax, : , 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold on;
title('x = 0.25 m');
grid on;
xlabel('Time (s)');
ylabel('Displacement (m)');


nexttile;
% C = 0.25
plot( t, U(0.5/deltax, :, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
plot( t, U(0.5/deltax, :, 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
plot( t, U(0.5/deltax, :, 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
plot( t, U(0.5/deltax, :, 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold on;
title('x = 0.5 m');
grid on;
xlabel('Time (s)');
ylabel('Displacement (m)');


nexttile;
% C = 0.25
plot( t, U(0.75/deltax, :, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
plot( t, U(0.75/deltax, :, 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
plot( t, U(0.75/deltax, :, 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
plot( t, U(0.75/deltax, :, 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold off;
title('x = 0.75 m');
grid on;
xlabel('Time (s)');
ylabel('Displacement (m)');


% Add legend to plot
h = [p1 p2 p3 p4];
lg = legend(h, {'C = 0.25','C = 0.5','C = 0.75','C = 1'},'Orientation', 'Horizontal');
lg.Layout.Tile = 'south';


%% Wave Ampltidue WTR x-coordinate
figure('Name','x-Dimension Plots','NumberTitle','off');
xp = tiledlayout(2,2); % 2 rows, 2 columns
xp.TileSpacing = 'compact';
xp.Padding = 'compact';
title(xp, 'Oscillating Displacement WRT x-Coordinate for Various C and t')
nexttile;
% C = 0.25
plot(x, U(:,1, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
plot( x, U(:,1, 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
plot( x, U(:, 1, 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
plot( x, U(:, 1, 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold on;
title('t = 0 s');
grid on;
xlabel('x-Coordinate (m)');
ylabel('Displacement (m)');
ylim([-1 1])

nexttile;
% C = 0.25
p1 = plot( x, U(:,1/deltat, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
p2 = plot( x, U(:,1/deltat, 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
p3 = plot( x, U(:, 1/deltat, 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
p4 = plot( x, U(:, 1/deltat, 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold on;
title('t = 1 s');
grid on;
xlabel('x-Coordinate (m)');
ylabel('Displacement (m)');


nexttile;
% C = 0.25
plot( x, U(:,2/deltat, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
plot( x, U(:,2/deltat, 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
plot( x, U(:, 2/deltat, 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
plot( x, U(:, 2/deltat, 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold on;
title('t = 2 s');
grid on;
xlabel('x-Coordinate (m)');
ylabel('Displacement (m)');


nexttile;
% C = 0.25
plot( x, U(:,3/deltat, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
plot( x, U(:,3/deltat, 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
plot( x, U(:, 3/deltat, 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
plot( x, U(:, 3/deltat, 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold off;
title('t = 3 s');
grid on;
xlabel('x-Coordinate (m)');
ylabel('Displacement (m)');
% Add legend to plot
h = [p1 p2 p3 p4];
lg = legend(h, {'C = 0.25','C = 0.5','C = 0.75','C = 1'},'Orientation', 'Horizontal');
lg.Layout.Tile = 'south';


%% Position dependant boundary condition 
figure('Name','x-Dimension BC','NumberTitle','off');
% C = 0.25
plot(x, U(:,1, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
plot( x, U(:,1, 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
plot( x, U(:, 1, 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
plot( x, U(:, 1, 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold on;
title('Position Dependant Boundary Condition for x-Domain');
grid on;
xlabel('x-Coordinate (m)');
ylabel('Displacement (m)');
ylim([-1 1])
legend('C = 0.25','C = 0.5','C = 0.75','C = 1', 'Location', 'southoutside', ...
    'Orientation', 'horizontal')

%% Time Dependant BC
figure('Name','Time Dep. BC.','NumberTitle','off');
% C = 0.25
plot( t, U(1,:, 0.25/Cstep), 'o', 'MarkerSize', 4, 'color', 'g', ...
        'LineStyle','--');
hold on;
% C = 0.50
plot( t, U(1,:, 0.50/Cstep), '*', 'MarkerSize', 4, 'color', 'b', ...
        'LineStyle','--');
hold on;
% C = 0.75
plot( t, U(1, :, 0.75/Cstep), '.', 'MarkerSize', 4, 'color', 'c', ...
        'LineStyle','--');
hold on;
% c = 1
plot( t, U(1, :, 1/Cstep), '+', 'MarkerSize', 4, 'color', 'r', ...
        'LineStyle','--');
hold on;
title('Time Dependant Boundary Condition for t0.01-Domain');
legend('C = 0.25','C = 0.5','C = 0.75','C = 1', 'Location', 'southoutside', ...
    'Orientation', 'horizontal');
grid on;
xlabel('Time (s)');
ylabel('Displacement (m)');
ylim([-1 1])