%% Incremental POD - Settings
tic
clc; clear; close all;
addpath('scripts/');
addpath('data/Cases/');

Only_liquid_phase = true;

%% LOAD DATA
[S, filename] = loadData(fullfile('Cases/','L8_G12.mat'));
frames = size(S.all_u_matrix_liquid, 3);

%% Vectorisation


[rows_liquid, columns_liquid, ~] = size(S.all_u_matrix_liquid);
Spatial_points_liquid = rows_liquid * columns_liquid;

if Only_liquid_phase
  
  U2_Vectorised = zeros(Spatial_points_liquid, frames);
  V2_Vectorised = zeros(Spatial_points_liquid, frames);
  parfor frame = 1:frames
    
    U2_Vectorised(:, frame) = reshape(S.all_u_matrix_liquid(:, :, frame), [], 1);
    V2_Vectorised(:, frame) = reshape(S.all_v_matrix_liquid(:, :, frame), [], 1);
    
  end
  
  Snapshot_matrix = [U2_Vectorised; V2_Vectorised];
  
else
  
  [rows_air, columns_air, ~] = size(S.all_u_matrix_air);
  Spatial_points_air = rows_air * columns_air;
  
  U2_Vectorised = zeros(Spatial_points_liquid, frames);
  V2_Vectorised = zeros(Spatial_points_liquid, frames);
  U1_Vectorised = zeros(Spatial_points_air, frames);
  V1_Vectorised = zeros(Spatial_points_air, frames);
  
  parfor frame = 1:frames
    
    U2_Vectorised(:, frame) = reshape(S.all_u_matrix_liquid(:, :, frame), [], 1);
    V2_Vectorised(:, frame) = reshape(S.all_v_matrix_liquid(:, :, frame), [], 1);
    U1_Vectorised(:, frame) = reshape(S.all_u_matrix_air(:, :, frame), [], 1);
    V1_Vectorised(:, frame) = reshape(S.all_v_matrix_air(:, :, frame), [], 1);
    
  end
  
  Snapshot_matrix = [U2_Vectorised; V2_Vectorised; U1_Vectorised; V1_Vectorised];
  
end


Snapshot_matrix(isnan(Snapshot_matrix)) = 0;

%% Normalisation

mean_matrix = mean(Snapshot_matrix, 2);

snapshot_fluctuations = Snapshot_matrix - mean_matrix;

%% Segment data

block_size = 350; % loading 100 frames at a time
number_of_blocks = ceil(frames / block_size); % calcualtes how many blocks we will process in total

snapshot_blocks = cell(number_of_blocks, 1);

parfor block = 1:number_of_blocks
  
  start_frame = ((block - 1) * block_size) + 1;
  end_frame = min(block * block_size, frames);
  snapshot_blocks{block} = snapshot_fluctuations(:, start_frame:end_frame);
  
end

%% Incremental POD

Temporal_Covariance = zeros(frames, frames);


for block = 1:number_of_blocks
  
  start_block = ((block - 1) * block_size) + 1;
  end_block = min(block * block_size, frames);
  current_block = snapshot_blocks{block};
  
  for upper_block = block:number_of_blocks
    
    start_upper_block = ((upper_block - 1) * block_size) + 1;
    end_upper_block = min(upper_block * block_size, frames);
    current_upper_block = snapshot_blocks{upper_block};
    
    computed_block = current_block' * current_upper_block;
    Temporal_Covariance(start_block:end_block, start_upper_block:end_upper_block) = computed_block;
    
    if block ~= upper_block
      
      Temporal_Covariance(start_upper_block:end_upper_block, start_block:end_block) = computed_block';
      
    end
  end
end

Temporal_Covariance = Temporal_Covariance / frames;


%% Eigenvalue Decompositon

[eigenvectors_matrix, eigenvalues_matrix] = eig(Temporal_Covariance);
eigenvalues_diagonal = diag(eigenvalues_matrix);
[eigenvalues, sort_index] = sort(eigenvalues_diagonal, 'descend');

sorted_eigenvalues = sort(eigenvalues, 'descend');

total_energy = sum(eigenvalues);
cumulative_energy = cumsum(eigenvalues) / total_energy;
thresholds = [0.90, 0.95, 0.99];
modes_to_retain = zeros(size(thresholds));

parfor threshold = 1:length(thresholds)
  
  modes_to_retain(threshold) = find(cumulative_energy >= thresholds(threshold), 1);
  
end

%% Calcuate spaital modes
number_of_modes_at_99 = find(cumulative_energy >= 0.99, 1);
Spatial_modes = snapshot_fluctuations * eigenvectors_matrix(:, sort_index(1:number_of_modes_at_99));

for mode = 1:number_of_modes_at_99
  Spatial_modes(:,mode) = Spatial_modes(:,mode) / norm(Spatial_modes(:,mode));
  
end

%% Saving results

results_directory = 'Results/POD_data';
if ~exist(results_directory, 'dir'), mkdir(results_directory); end
[~, base_name, ~] = fileparts(filename);
save_path = fullfile(results_directory, ['POD_Results_', base_name, '.mat']);
save(save_path, 'mean_matrix', 'eigenvectors_matrix','eigenvalues', 'Spatial_modes','sort_index','cumulative_energy','rows_liquid', 'columns_liquid', 'Only_liquid_phase', 'filename','number_of_modes_at_99', '-v7.3');

%% Plotting
%--TODO: Add x and y labels as well as a legend

plot(1:length(cumulative_energy), cumulative_energy(:)' * 100, 'b-', 'Linewidth', 2)
hold on;

for threshold = 1:length(thresholds)
  
  yline(thresholds(threshold) * 100, '--');
  xline(modes_to_retain(threshold), ':');
  plot(modes_to_retain(threshold), thresholds(threshold) * 100, 'o');
  
end

%% Plot time coefficinets for the first 10 modes
%--TODO: Plot time coefficinets of the first 10 modes,
%--TODO: Identify the pairs (look for sine/cose sine) that are out of phase
%--TODO: If they look too messy run the hilbert() command to identify pairs

figure('Name', 'Time Coefficients for the first 10 modes')
frames_per_second = 1000;
time_frame = ((1:frames)/1000);

for mode =  1:10
  mode_index = sort_index(mode);
  subplot(5, 2, mode)
  plot(time_frame,eigenvectors_matrix(:,mode_index));
  title('Mode', mode);
  xlabel('Frame')
  ylabel('Mode Coefficient')
  xlim([0, 1]);
  yticks([-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08])
end
hold off;
%% identify pairs via phase plots:

% Hilbert Transforms:

Hilbert_mode_a = hilbert(eigenvectors_matrix(:, sort_index(mode_a)));
Hilbert_mode_b = hilbert(eigenvectors_matrix(:, sort_index(mode_b)));

phase_diffrence = rad2deg(angle(Hilbert_mode_a ./ Hilbert_mode_b));
average_phase_shift = median(phase_diffrence);
disp(average_phase_shift);

% Plotting pairs
figure('Name', 'Phase plots comparisons')

mode_a = 1;
mode_b = 2;

plot(time_frame, eigenvectors_matrix(:, sort_index(mode_a)), 'DisplayName', ['Mode', num2str(mode_a)]);
hold on;
plot(time_frame, eigenvectors_matrix(:, sort_index(mode_b)), 'DisplayName', ['Mode', num2str(mode_b)]);
xlabel('Frame')
ylabel('Mode Coefficient')
xlim([0, 0.5]);
yticks([-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08])
title(sprintf('Mode %d Vs Mode %d with a phase shift of %.3f for case: %s', mode, mode + 1, average_phase_shift, filename));
legend('FontSize', 9)

toc

