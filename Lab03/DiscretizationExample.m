%% Example from Lab 3
%Clear commands and workspace
clc; 
clear all;
close all; 

% Variable declaration
Nx = 10; % Number of grid points in X
Ny = 10; % Number of grid points in Y

x = 1:1:Nx; % x-grid, delta-x
y = 1:1:Ny; % y-grid, delta-y
[X, Y] = meshgrid(x, y); % Create a 2D grid

plot(X, Y, 'ko', 'MarkerFaceColor', 'r'); % Plot with trd points