%% Incremental POD

clc, clear, close all,


%% LOAD DATA

[S, filename] = loadData('L8_G7.mat');

[X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);
% 1's are the air, 2's are the liquid

frames = size(S.all_u_matrix_liquid, 3);

%% Vectorisation
% POD Requires 2D matrices, not 3D blocks.
% So Every frame of the flow field must be flatterned from a 2D grid into a single, 1D column vector

% obtain the rows and columns of the liquid to setup matrices
[rows_liquid, columns_liquid, ~] = size(U2); 
Spatial_points = rows_liquid * columns_liquid;

% initalise matrices to store vectorised velocity components
U2_Vectorised = zeros(Spatial_points, frames);
V2_Vectorised = zeros(Spatial_points, frames);






