%% Incremental POD

clc, clear, close all,kk


%% LOAD DATA

[S, filename] = loadData('L8_G7.mat');

[X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);
% 1's are the air, 2's are the liquid

%% Vectorisation
% POD Requires 2D matrices, not 3D blocks.
% So Every frame of the flow field must be flatterned from a 2D grid into a single, 1D column vector

[rows_liquid, columns_liquid, ~] = size(U2) 




