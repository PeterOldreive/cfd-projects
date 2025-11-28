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
omega          = 1.95;     % relaxation factor (ω = 1 → pure Gauss–Seidel)
tol            = 1;    % convergence tolerance 
maxIterations = 5000; % Maximum number of iterations before forcing 
                       % convergance of the algorithm
guess = 1000; % First guess of solution vector values 

%% Define function for the element wise Gauss-Siedel Formula
function xik1 = element_based(xj_1, xj, aij, bi, aii, index)
    % This function solves the element based equation for the Gauss-Seidel
    % method fof element xi(k+1), using xj(k+1), xj, aij, bi
    % Find sums of aijxj(k+1) and aijxj(k)
    
    sumxj1 = 0;
    % Sum aij, xj(k+1)
    for i = 1:(index - 1)
        sumxj1 = sumxj1 + aij(i)*xj(i);
    end 

    sumxj = 0; 
    % Sum aij, xj
    for i = (index+1):length(xj_1)
        sumxj = sumxj + aij(i)*xj_1(i);
    end 

    % Apply the element based formula 
    xik1 = (bi - sumxj1 - sumxj) / aii;
end

%% Import .csv Files from External Script 

% % Import matrx of coefficents
% coeffMatrix = readmatrix('coeff_matrix.csv');
% % Import vector of knowns
% knownVector = readmatrix('knowns.csv');
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
for iter = 1:maxIterations % Iterate until maximum number has been reached
    xj_1 = xj; 
    disp(iter)
    for i = 1:coeffDim(1) % Loop for number of rows in coeff matrix
        aii = coeffMatrix(i, i); % Coefficent aii lies along the main 
                                 % diagonal of the matrix of coefficents 
        aij = coeffMatrix(i, :); % all coefficents for the row i in the 
                                 % matrix of coefficents
        bi = knownVector(i);   % Pull appropriate known value from vector 
        
        % Call function to impliment elementwise formula 
        %xj(i) = element_based(xj_1, xj, aij, bi, aii, i); 
        sumxj1 = 0;
        % Sum aij, xj(k+1)
        for j = 1:(i - 1)
            sumxj1 = sumxj1 + aij(j)*xj(j);
        end 

        sumxj = 0; 
        % Sum aij, xj
        for j = (i+1):length(xj_1)
            sumxj = sumxj + aij(j)*xj_1(j);
        end 

        % Apply the element based formula 
        xj(i) = (bi - sumxj1 - sumxj) / aii;
        % SOR update
        xj(i) = (1-omega)*xj_1(i) + omega * xj(i);
            
    end
    % Check for convergance 
    if max(abs(xj - xj_1)) < tol
        % If the maximum absulte difference between all new and old
        % vaules is less than the tolarance, break out of the loop. 
        fprintf('Converged after %d iterations.\n', iter);
        break; % Exit the loop if converged
    end
end

% %% Output results 
% % Save solution vector to a .csv file
% writematrix(xj, 'solution_vector.csv');

return; 