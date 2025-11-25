
clear; clc; close all;

%% Geometry: 4 x 4 nodes 

Nx = 4;              % nodes in x
Ny = 4;              % nodes in y
N_total = Nx * Ny;   % total number of nodes

% Allocate coefficient arrays
N  = zeros(N_total,1);   % diagonal coefficients a_P
NE = zeros(N_total,1);   % east  neighbour
NW = zeros(N_total,1);   % west  neighbour
NN = zeros(N_total,1);   % north neighbour
NS = zeros(N_total,1);   % south neighbour
B  = zeros(N_total,1);   % source term (right-hand side)

% Boundary temperatures (Dirichlet)
T_left  = 100;   % left wall
T_right = 50;    % right wall
T_top   = 50;    % top wall
T_bottom= 50;    % bottom wall

%% Set up the coefficients using a finite-volume style loop

for i = 1:N_total
    
    row    = floor((i-1)/Nx) + 1;   % j = 1..Ny (y index)
    column = mod((i-1),Nx) + 1;     % k = 1..Nx (x index)
    
    % ---- Left boundary (T = 100) ----%
    if column == 1
        N(i) = 1;          % a_P
        B(i) = T_left;     % RHS
        % NE,NW,NN,NS remain 0
        
    % ---- Right boundary (T = 50) ----%
    elseif column == Nx
        N(i) = 1;
        B(i) = T_right;
        
    % ---- Bottom boundary (T = 50) ----%
    elseif row == 1
        N(i) = 1;
        B(i) = T_bottom;
        
    % ---- Top boundary (T = 50) ----%
    elseif row == Ny
        N(i) = 1;
        B(i) = T_top;
        
    % ---- Interior nodes: Laplace stencil ----%
    else
        N(i)  = 4;     % a_P
        NE(i) = -1;
        NW(i) = -1;
        NN(i) = -1;
        NS(i) = -1;
        B(i)  = 0;
    end
end

%% Gauss–Seidel with SOR

omega          = 2;     % relaxation factor (ω = 1 → pure Gauss–Seidel)
tol            = 1e-6;    % convergence tolerance
max_iterations = 5000;

% initial guess (say, uniform 50)
T = 50 * ones(N_total,1);

for iter = 1:max_iterations
    
    T_old = T;
    
    for i = 1:N_total
        
        % --- Compute Gauss–Seidel value T_GS(i)
        if i == 1
            % bottom-left corner (has east and north neighbour only)
            T_GS = ( B(i) - NE(i)*T(i+1)   - NN(i)*T(i+Nx) ) / N(i);
            
        elseif i > 1 && i <= Nx
            % bottom edge (excluding i=1)
            T_GS = ( B(i) - NE(i)*T(i+1)   - NW(i)*T(i-1) ...
                           - NN(i)*T(i+Nx) ) / N(i);
            
        elseif i > N_total-Nx && i < N_total
            % top edge (excluding top-right corner)
            T_GS = ( B(i) - NE(i)*T(i+1)   - NW(i)*T(i-1) ...
                           - NS(i)*T(i-Nx) ) / N(i);
            
        elseif i == N_total
            % top-right corner (has west and south neighbour only)
            T_GS = ( B(i) - NW(i)*T(i-1)   - NS(i)*T(i-Nx) ) / N(i);
            
        else
            % true interior nodes
            T_GS = ( B(i) - NE(i)*T(i+1)   - NW(i)*T(i-1) ...
                           - NN(i)*T(i+Nx) - NS(i)*T(i-Nx) ) / N(i);
        end
        
        % --- SOR update ---
        T(i) = (1-omega)*T(i) + omega * T_GS;
    end
    
    % convergence check
    if max(abs(T - T_old)) < tol
        fprintf('Converged in %d iterations, max change = %.3e\n', ...
                iter, max(abs(T-T_old)));
        break;
    end
end

%%   Post-processing: reshape and inspect interior unknowns

T_field = reshape(T, [Nx, Ny])';  
disp('Temperature field (rows = y, cols = x):');
disp(T_field);

% Interior unknowns 
T11 = T_field(2,2);
T21 = T_field(2,3);
T12 = T_field(3,2);
T22 = T_field(3,3);

fprintf('T11 = %.5f\n', T11);
fprintf('T21 = %.5f\n', T21);
fprintf('T12 = %.5f\n', T12);
fprintf('T22 = %.5f\n', T22);

% Optional contour plot (For practice only)
x = 0:Nx-1;
y = 0:Ny-1;
figure;
contourf(x, y, T_field, 10, 'k');
colorbar; colormap('jet');
xlabel('x index');
ylabel('y index');
title(sprintf('Steady 2D conduction, \\omega = %.2f', omega));
axis equal tight;
