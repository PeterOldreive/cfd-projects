%% Problem with 2D velocity contour
clc; clear all; close all;

% variable declaration
Nx = 10;           % Number of grid points in x-direction
Ny = 10;           % Number of grid points in y-direction

x = 1:1:Nx;          % x-grid, delx=1
y = 1:1:Ny;          % y-grid, dely=1

[X, Y] = meshgrid(x, y);   % 2D grid coordinates
% 
plot(X, Y, 'ko', 'MarkerFaceColor', 'r'); % red points
grid on;

U = randn(length(y), length(x));  % match X,Y      

% % % ---- Plot the 2D field ----
figure;
contourf(X, Y, U, 20);     % Filled contour plot with 0levels
colorbar;                  % Show color scale
title('2D Velocity Field');
xlabel('x-axis');
ylabel('y-axis');
grid on;
%axis equal tight;          % Equal scaling and tight axes

hold on;

% Show grid points as red dots
plot(X, Y, 'ro', 'MarkerSize', 6, 'LineWidth', 1);

% Show grid lines
% 
for i = 1:Ny
      plot(x, Y(i,:), 'k-', 'LineWidth', 0.5); % horizontal lines
end

for j = 1:Nx
   plot(X(:,j), y, 'k-', 'LineWidth', 0.5); % vertical lines
end

% hold off;

%% Find diagonal and diagonal matrix 
 n = size(U,1);              % number of rows (since U is square)

% % % diagonal values
diag_vals = [];             % empty vector

for i = 1:n 
     diag_vals(i) = U(i,i);  % row = column
end

% % % diagonal matrix

D = zeros(n,n);             % start with zeros
for i = 1:n
      D(i,i) = diag_vals(i);  % put values on diagonal
end

disp('Diagonal values:');
disp(diag_vals)

% % disp('Diagonal matrix:');
% % disp(D)

%% Plot diagonal values vs x-axis 
x = 1:length(diag_vals);           % x-axis(1 to n)

% %y = 1:length(diag_vals); 

% figure;
plot(x, diag_vals, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
% % %plot(diag_vals, y, 'ro--', 'LineWidth', 1.5, 'MarkerSize', 6);

% grid on;
% % title('Diagonal Elements vs X-axis');
% xlabel('x axis');
% ylabel('Diagonal velocity');

%% plot 5th row w. r. to x axis of your matris
row5 = U(5, :);
disp('5th row of U:');
disp(row5);

x = 1:size(U,2);           % x-axis for columns
figure;
plot(x, row5, 'r-o', 'LineWidth', 1.5);
xlabel('Column index (x)');
ylabel('U(5, :) values');
grid on;

%%
% Excercise: plot 6th column of a matrix U with respect to 
% y-axis.

%% control --> if else, for loop, nested loop

% Find positive and negative

clc; clear all; close all;

Nx = 10; Ny = 10;

x = 1:Nx;
y = 1:Ny;

U = randn(length(y), length(x)); 

disp('Matrix U:')
disp(U)

% ---- If-else condition on each element/node
for i = 1:Ny
    for j = 1:Nx
        if U(i,j) > 0
            fprintf('U(%d,%d) = %.2f is Positive\n', i, j, U(i,j));
        else
            fprintf('U(%d,%d) = %.2f is Negative\n', i, j, U(i,j));
        end
    end
end

% if elseif else:
% for i = 1:Ny
%     for j = 1:Nx
%         if U(i,j) > 0
%             fprintf('U(%d,%d) = %.2f is Positive\n', i, j, U(i,j));
%         elseif U(i,j) == 0
%             fprintf('U(%d,%d) = %.2f is Zero\n', i, j, U(i,j));
%         else
%             fprintf('U(%d,%d) = %.2f is Negative\n', i, j, U(i,j));
%         end
%     end
% end

%% Exercise:
% 1. Data types (double, int, float, string, ...)
% 2. Arrays (list, tuple, dictionary)
% 3. Control-(For loop, while loop, if else conditions)
% 4. Matrix - indexing, slicing ...

%%%%%%%%%%%%%%%%%%%%%%%%%% By Khalid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%