clear; clc; close all

% Load data
S = loadData('L8_G3.mat');  % Pass the filename directly as a string
frame = 10;

%Get all data at once
[X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);

%% Create masks
[Z1_masked, Z2_masked] = createMasks(X1, Y1, X2, Y2, Z1, Z2, X3, Y3);

%% Plot Velocity Vectors with heatmap
plotVectorMap(X1, Y1, U1, V1, X2, Y2, U2, V2, X3, Y3, Z1_masked, Z2_masked);

%% Animate Velocity Vectors
frames = 100;
animateVelocityVectors(S, frames);
