clc; clear; close all;

[S, filename] = loadData('L8_G9.mat'); 
frames = size(S.smoothed_film_height_matrix_out, 2) - 1;

heights = [];

for frame = 1:frames
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Y3] = getData(S, frame);
    heights = [heights; Y3(:)];
end


% stats
h_max = max(heights);
h_min = min(heights);
h_mean = mean(heights);

fprintf('Film Height Statistics (All %d Frames):\n', frames);
fprintf('Total Data Points: %d\n', length(heights));
fprintf('Max Height: %.2f mm\n', h_max);
fprintf('Min Height: %.2f mm\n', h_min);
fprintf('Mean Height: %.2f mm\n', h_mean);

bin_edges = [h_min, 1.0, 1.5, 2.0, 2.5, 3.5, h_max];
n_bins = length(bin_edges) - 1;

% count data points in each bin
bin_counts = zeros(n_bins, 1);
fprintf('\nData Points in Each Bin:\n');

for bin = 1:n_bins
    bin_counts(bin) = sum(heights >= bin_edges(bin) & heights < bin_edges(bin+1));
    fprintf('Bin %d (%.2f - %.2f mm): %d data points\n', bin, bin_edges(bin), bin_edges(bin+1), bin_counts(bin));
end

% create storage for vertical velocity data in each bin

bin_data = repmat(struct( ...
    "U1", [], ...
    "V1", [], ...    
    "U2", [], ...
    "V2", [] ...
), n_bins, 1);

Y_profile = []; % store vertical positions

% populate bins with vertical profiles based on film height

tic
for frame = 1:1000

    fprintf('Processing frame %d / %d\r', frame, frames);
    [X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);

    [Ny1, Nx1] = size(U1);
    X_air_columns = X1(1, :);

    if isempty(Y_profile)
        Y_profile = Y1(:, 1); % storing once since y is consistent across x 
    end
    
    % air phase
    for col = 1:Nx1
        x_pos = X_air_columns(col);
    end

    % local film height at this frame
    local_film_height = interp1(X3, Y3, x_pos, 'linear', 'extrap');
    if isnan(local_film_height)
        continue; % skip if no valid film height
    end
    % determine bin
    bin_index = find(local_film_height >= bin_edges(1:end-1) & local_film_height < bin_edges(2:end), 1);

    if isempty(bin_index) && local_film_height == bin_edges(end)
        bin_index = n_bins; % assign to last bin if equal to max edge
    elseif isempty(bin_index)
        continue; % skip if no bin found
    end

    % extract vertical profile for this column
    U1_column = U1(:, col);
    V1_column = V1(:, col);

    % append as new column in chosen bin
    bin_data(bin_index).U1 = [bin_data(bin_index).U1, U1_column];
    bin_data(bin_index).V1 = [bin_data(bin_index).V1, V1_column];

    % liquid phase
    [Ny2, Nx2] = size(U2);
    X_liquid_columns = X2(1, :);

    for col = 1:Nx2
        x_pos = X_liquid_columns(col);

        % local film height at this frame
        local_film_height = interp1(X3, Y3, x_pos, 'linear', 'extrap');
        if isnan(local_film_height)
            continue; % skip if no valid film height
        end
        % determine bin
        bin_index = find(local_film_height >= bin_edges(1:end-1) & local_film_height < bin_edges(2:end), 1);

        if isempty(bin_index) && local_film_height == bin_edges(end)
            bin_index = n_bins; % assign to last bin if equal to max edge
        elseif isempty(bin_index)
            continue; % skip if no bin found
        end

        % extract vertical profile for this column
        U2_column = U2(:, col);
        V2_column = V2(:, col);

        % append as new column in chosen bin
        bin_data(bin_index).U2 = [bin_data(bin_index).U2, U2_column];
        bin_data(bin_index).V2 = [bin_data(bin_index).V2, V2_column];
    end 
end 
toc
