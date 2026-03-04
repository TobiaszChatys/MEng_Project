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
% U2 and V2_Vectorised are 4950 by 2499 meaning once we create the snapshot matrix it should be 9900 by 2499

% Create the snapchot matrix
Snapshot_matrix = [U2_Vectorised; V2_Vectorised];
% Snapshot matrix is 9900 by 2499

% removes all nans from the matrix
Snapshot_matrix(isnan(Snapshot_matrix)) = 0;

%% Normalisation
% A flow field can be split into a steady, time average component u dash, and a fluctuating dynamic
% component u', POD is always perfomed on the fluctuations.
% By subtracting the time-averaged mean from every single frame, the POD will isolate the moving,
% turbulent structures rather than just showing static background noise. centering the data, ensuring the
% resulting POD modes describe how the flow changes or fluctuates over time.

% taking the mean across every column, mean(snapshot_matrix, 1) would take the mean at every row
mean = mean(Snapshot_matrix, 2); % 9900 by 1 matrix needs to be replicated

% computing the fluctuations u', altough mean is a 9900 by 1 matrix and snapshot_matrix is 9900 by 2499
% MATLAB will notice that it only has one column an preten it has been repeated to match the snapshot_matrix
snapshot_fluctuations = Snapshot_matrix - mean;

%% Segment data
% Because the fluctuations matrix X is huge, loading it all into the memory for matrix multiplication can 
% slow down the process.
% so we can process the matrix by loading block by block

block_size = 100; % loading 100 frames at a time
number_of_blocks = ceil(frames / block_size); % calcualtes how many blocks we will process in total
% ceil rounds up to make sure there are no left over frames 

% A cell array makes sure that if the last block is smaller than the others, the system wont crash. A cell
% can hold diffrent sizes of matries whereas normal matrices must be perfectly rectangular.
Snapshot_blocks = cell(number_of_blocks, 1);

% We need to perform some logic to tell MATLAB where to start and finish each block
% it will always be (n * 100) + 1
parfor block = 1:number_of_blocks
    start_frame = ((block - 1) * block_size) + 1;
    end_frame = min(block * block_size, frames);
    snapshot_blocks{block} = snapshot_fluctuations(:, start_frame:end_frame);
end









