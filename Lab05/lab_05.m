clc; clear all;

% Load Data coming from Ax=b 
load('Data_A_x_b_and_T.mat')

%% Check the residual error (to verify solution accuracy)
residual = norm(A*x - b);

% Display results
fprintf('Residual error: %.2e\n', residual);

%% Visualize of 2D for matrix T
figure;
imagesc(T_plot);
colormap('jet');
colorbar;
title('Temperature distribution');
xlabel('x-axis');
ylabel('y-axis');

%% visualization of 3D for matrix T
figure;
surf(T_plot, 'EdgeColor','k');
xlabel('x-axis'); 
ylabel('y-axis');
zlabel('Temperature distribution');
colormap('jet');
colorbar; 
view(40,30);

%% Plot diagonal values vs x-axis 
diag_vals_T = diag(T_plot);     % main diagonal of T
n = length(diag_vals_T);   % number of diagonal points

%x = 1:n;  % x-axis index (1 to n)
y = 1:n;  % y-axis index (1 to n)

figure;
%plot(x, diag_vals_T, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(diag_vals_T, y, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
title('Diagonal elements vs x-axis');
xlabel('x axis');
ylabel('Diagonal velocity');

%% plot 25th row with respect to x axis of your matrix

row25 = T_plot(25, :); % extract the 25th row

x = 1:size(T_plot, 2); % x-axis (column indices)
figure;
plot(x, row25, 'r-o', 'LineWidth', 1.5);
xlabel('x-axis');
ylabel('T(25, :) values');
title('25th Row of T vs. x-axis');
grid on;

% If you want to see the data: 
disp('25th row of T_plot:');
disp(row25);

%% plot 25th column with respect to y axis of your matrix

col25 = T_plot(:, 25); % extract the 25th row

y = 1:size(T_plot, 2); % x-axis (column indices)
figure;
plot(col25, y, 'r-o', 'LineWidth', 1.5);
xlabel('T(:, 25) values');
ylabel('y-axis');
title('25th colmun T vs. y-axis');
grid on;

%% Contour plot

[Nx, Ny] = size(T_plot);
x = linspace(0, 1, Nx);   % adjust your step size
y = linspace(0, 1, Ny);
[X, Y] = meshgrid(x, y);

figure;
contourf(X, Y, T_plot, 50, 'LineColor', 'none');  % 50 = number of contour levels
colorbar;
colormap('jet');
title('Contour plot of Temperature');
xlabel('x');
ylabel('y');
axis equal
grid on

%%
[Ux, Vy] = gradient(T_plot);
figure;
contourf(X, Y, T_plot, 50, 'LineColor', 'none');
hold on;
quiver(X, Y, Ux, Vy, 'w', 'LineWidth', 2);
colorbar;
colormap('jet');
xlabel('x-axis');
ylabel('y-axis');
title('Gradient Field of Temperature');
axis equal;
grid on;

%%

% --- Compute gradient field (acts as velocity field)
[Ux, Vy] = gradient(T_plot);

% Compute magnitude for background and streamline coloring
mag = sqrt(Ux.^2 + Vy.^2);

% --- Create filled contour background of scalar field
figure('Color','w');
contourf(X, Y, T_plot, 60, 'LineColor', 'none');
hold on;

% --- Overlay smooth streamlines
% Define starting points 
[startX, startY] = meshgrid( ...
    linspace(min(X(:)), max(X(:)), 25), ...   % number of seed points along X
    linspace(min(Y(:)), max(Y(:)), 25));      % number of seed points along Y

% Draw streamlines
streamObj = streamline(X, Y, Ux, Vy, startX, startY);
set(streamObj, ...
    'Color', 'w', ...        % streamline color (white stands out on jet)
    'LineWidth', 1.5);       % thickness for visibility


colorbar;
colormap('jet');             
xlabel('x-axis');
ylabel('y-axis');
title('Streamlines of Gradient Field (Temperature)');
axis equal tight;
box on;
grid on;

%% Histogram plot
figure('Color','w');
histogram(T_plot(:), 20, 'FaceColor','b', 'EdgeColor','k');
xlabel('T values');
ylabel('Frequency');
title('Temperature Distribution Histogram');
grid on;

%%%%%%%%%%%%%% By Khalid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















