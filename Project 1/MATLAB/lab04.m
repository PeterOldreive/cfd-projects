clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Boundary values (Dirichlet BCs)

T_left   = 100;   % left wall
T_right  = 50;    % right wall
T_top    = 50;    % top wall
T_bottom = 50;    % bottom wall

% Matrix declaration
A = zeros(4,4);
b = zeros(4,1);

%% Equation for T11  -> row 1
A(1,1) = -4;   % T11 coefficient
A(1,2) =  1;   % T21 coefficient
A(1,3) =  1;   % T12 coefficient
A(1,4) =  0;   % T22 coefficient
b(1)   = -(T_left + T_bottom);   % -100 - 50 = -150

%%
% Equation for T21  -> row 2
A(2,1) =  1;   % T11
A(2,2) = -4;   % T21
A(2,3) =  0;   % T12
A(2,4) =  1;   % T22
b(2)   = -(T_right + T_bottom);  % -50 + 50 = -100
%%
% Equation for T12  -> row 3
A(3,1) =  1;   % T11
A(3,2) =  0;   % T21
A(3,3) = -4;   % T12
A(3,4) =  1;   % T22
b(3)   = -(T_left + T_top);      % -100 + 50 = -150
%%
% Equation for T22  -> row 4
A(4,1) =  0;   % T11
A(4,2) =  1;   % T21
A(4,3) =  1;   % T12
A(4,4) = -4;   % T22
b(4)   = -(T_right + T_top);     % -50 + 50 = -100

%% Solve equation
x = A \ b;

T11 = x(1);   % bottom-left interior
T21 = x(2);   % bottom-right interior
T12 = x(3);   % top-left interior
T22 = x(4);   % top-right interior

%% Reconstruct full field including the wall temperatures

corner_top_left     = T_left;    % 100
corner_top_right    = T_right;   % 50
corner_bottom_left  = T_left;    % 100
corner_bottom_right = T_right;   % 50

T_full = [ corner_top_left,   T_top,     T_top,     corner_top_right;   % top boundary row
           T_left,            T12,       T22,       T_right;            % interior top row
           T_left,            T11,       T21,       T_right;            % interior bottom row
           corner_bottom_left,T_bottom,  T_bottom,  corner_bottom_right]; % bottom boundary row

% For plotting with y increasing upward, flip vertically for imagesc later
T_plot = flipud(T_full);

%% 5. Post-process / output
fprintf('Solved interior temperatures:\n');
fprintf('  T11 (bottom-left interior)  = %.2f\n', T11);
fprintf('  T21 (bottom-right interior) = %.2f\n', T21);
fprintf('  T12 (top-left interior)     = %.2f\n', T12);
fprintf('  T22 (top-right interior)    = %.2f\n', T22);
%%
disp(' ');
disp('Full temperature field including boundaries (top row first):');
disp(T_full);

%%
figure;
surf(T_plot,'EdgeColor','k');
xlabel('x-index');
ylabel('y-index');
zlabel('Temperature');
title('Temperature field (surface view)');
colorbar;
view(45,45);
%%
figure;
imagesc(T_plot);
set(gca,'YDir','normal'); % make y increase upward visually
axis equal tight;
colorbar;
title('Temperature contour map');
xlabel('x-index');
ylabel('y-index');


