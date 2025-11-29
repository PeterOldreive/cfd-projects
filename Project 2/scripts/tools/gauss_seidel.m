% This script implements the Gauss-Seidel method for approximating
% the solution to a linear system of equations. 
%
% This script can be implemented into other scripts for solving systems
% using the Gauss-Seidel method 
%
% Author: Peter Oldreive
% Created: 2025-11-24
% Updated: 2025-11-24
%% Convergance specific parameters 
omega          = 1.97;     % relaxation factor (ω = 1 → pure Gauss–Seidel)
tol            = 0.001;    % convergence tolerance 
maxIterations = 20000; % Maximum number of iterations before forcing 
                       % convergance of the algorithm
guess = 800; % First guess of solution vector values 

%% Import .csv Files from External Script 

% % Import matrx of coefficents
% A = readmatrix('coeff_matrix.csv');
% % Import vector of knowns
% B = readmatrix('knowns.csv');
%% Take Matricies Directly from workspace
coeffMatrix = A; 
knownVector = B;
%% Checks
% Check that vector has the same mumber of elements as rows in the 
% coefficent matrix. If these are not compatible a warning will be given
if (length(coeffMatrix(:,1)) ~= length(knownVector)) | ...
    (length(coeffMatrix(1,:)) ~= length(knownVector))
    fprintf('Matrix and Vector not Compatble \n')
    % Write a 2x1 matrix of zeros to provide some output 
    writematrix(zeros(2,1), 'solution_vector.csv');
    return; % Pause running script 
end 

%% Initialize first Guess of Solution Vector
% Initialize the guess for temperature values (assume 300 K)
coeffDim = size(coeffMatrix); % Save dimensions of coefficent matrix
% Initialize solution vector with 300 K 
xj = guess * ones(coeffDim(1),1); 

%% Run iteritve Gauss-Seidel Method 
for iter = 1:maxIterations
    xj_1 = xj; 
    for i = 1:Ny
        index = (i-1)*Nx + 1;
        for j = 1:Nx
            term_E = 0;
            if(j < Nx)
                term_E = xj_1(index+1)*A(index, index + 1); 
            end
            term_W = 0;
            if(j > 1)
                term_W = xj(index-1)*A(index, index - 1); 
            end
            term_N = 0;
            if(i < Ny)
                term_N = xj_1(index+Nx)*A(index, index + Nx);
            end
            term_S = 0;
            if(i > 1)
                term_S = xj(index-Nx)*A(index, index - Nx);
            end
        
            xj(index) = (B(index) - term_N - term_S - term_E - term_W)/A(index,index);
            xj(index)= (omega)*xj(index) + (1-omega)*xj_1(index);
            index = index + 1; 
        end 
    end 
          % Check for convergance 
        if max(abs(xj - xj_1)) < tol
            % If the maximum absulte difference between all new and old
            % vaules is less than the tolarance, break out of the loop. 
            %fprintf('Converged after %d iterations.\n', iter);
            disp(['Converged in: ' num2str(iter) ' Iterations'])
            break; % Exit the loop if converged
        end 
end

% %% Output results 
% % Save solution vector to a .csv file
% writematrix(xj, 'solution_vector.csv');

return; 