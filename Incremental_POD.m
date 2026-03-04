%% Incremental POD

clc, clear, close all,


%% LOAD DATA

[S, filename] = loadData('L8_G7.mat'); 
frames = size(S.all_u_matrix_liquid, 3);

%% Vectorisation
% POD Requires 2D matrices, not 3D blocks.
% So Every frame of the flow field must be flatterned from a 2D grid into a single, 1D column vector

% obtain the rows and columns of the liquid to setup matrices
[rows_liquid, columns_liquid, ~] = size(S.all_u_matrix_liquid);
Spatial_points = rows_liquid * columns_liquid;

% initalise matrices to store vectorised velocity components
U2_Vectorised = zeros(Spatial_points, frames);
V2_Vectorised = zeros(Spatial_points, frames);


% reshape function: reshape(A, sz), reshapes A using size vector sz, to define the resultant
% the brackets tell matlab to automatically calculate how many rows are needed to fit all the data
parfor frame = 1:frames
    U2_Vectorised(:, frame) = reshape(S.all_u_matrix_liquid(:, :, frame), [], 1);
    V2_Vectorised(:, frame) = reshape(S.all_v_matrix_liquid(:, :, frame), [], 1);
end

size(V2_Vectorised)
size(U2_Vectorised)





