clc; clear; close all;

[S, filename] = loadData('L8_G9.mat'); 
frames = size(S.all_u_matrix_liquid, 3);

[
X1, Y1, U1, V1, Z1, ... % air phase 
X2, Y2, U2, V2, Z2,  ... % liquid phase 
X3, Y3  ... % film data
] = getData(S, 1);

all_film_heights = [];

for frame = 1:frames
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Y3] = getData(S, frame);
    all_film_heights = [all_film_heights; Y3(:)];
end

%% Compute film height statistics & define bins

film_max_height = max(all_film_heights);
film_min_height = min(all_film_heights);
film_mean_height = mean(all_film_heights);
film_mode_height = mode(all_film_heights);

fprintf('Film Height Statistics (All %d Frames):\n', frames);
fprintf('Total Data Points: %d\n', length(all_film_heights));
fprintf('Max Height: %.2f mm\n', film_max_height);
fprintf('Min Height: %.2f mm\n', film_min_height);
fprintf('Mean Height: %.2f mm\n', film_mean_height);
fprintf('Mode Height: %.2f mm\n', film_mode_height);

bin_edges = [film_min_height, 1.0, 1.5, 2.0, 2.5, 3.5, film_max_height];
number_of_bins = length(bin_edges) - 1;

% count data points in each bin
bin_counts = zeros(number_of_bins, 1);
fprintf('\nData Points in Each Bin:\n');

for bin = 1:number_of_bins
    bin_counts(bin) = sum(all_film_heights >= bin_edges(bin) & all_film_heights < bin_edges(bin+1));
    fprintf('Bin %d (%.2f - %.2f mm): %d data points\n', bin, bin_edges(bin), bin_edges(bin+1), bin_counts(bin));
end

%% create storage for vertical velocity data in each bin

bin_data = repmat(struct( ...
    "U1", [], ...
    "V1", [], ...    
    "U2", [], ...
    "V2", [] ...
), number_of_bins, 1);

Y_profile_air = Y1(:,1); % store air vertical positions
Y_profile_liquid = Y2(:,1); % store liquid vertical positions

[number_of_yair, number_of_xair] = size(U1); % getting dimensions
[number_of_yliquid, number_of_xliquid] = size(U2);

%% populate bins with vertical profiles based on film height

tic

for frame = 1:frames

    fprintf('Processing frame %d / %d\r', frame, frames);
    [X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);

    X_air_columns = X1(1, :);
    X_liquid_columns = X2(1, :);

    % air phase
    for column = 1:number_of_xair
        
        x_pos = X_air_columns(column);
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
        U1_column = U1(:, column);
        V1_column = V1(:, column);

        % append as new column in chosen bin
        bin_data(bin_index).U1 = [bin_data(bin_index).U1, U1_column];
        bin_data(bin_index).V1 = [bin_data(bin_index).V1, V1_column];
    end

    % liquid phase
    for column = 1:number_of_xliquid
        x_pos = X_liquid_columns(column);

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
        U2_column = U2(:, column);
        V2_column = V2(:, column);

        % append as new column in chosen bin
        bin_data(bin_index).U2 = [bin_data(bin_index).U2, U2_column];
        bin_data(bin_index).V2 = [bin_data(bin_index).V2, V2_column];
    end
end

toc

%% compute conditonal mean profiles for each bin

conditional_means = repmat(struct( ...
    "U1_mean", [], ...
    "V1_mean", [], ...
    "U2_mean", [], ...
    "V2_mean", [] ...
), number_of_bins, 1);

for bin = 1:number_of_bins
    if ~isempty(bin_data(bin).U1)
        conditional_means(bin).U1_mean = mean(bin_data(bin).U1, 2, 'omitnan');
        conditional_means(bin).V1_mean = mean(bin_data(bin).V1, 2, 'omitnan');
    end

    if ~isempty(bin_data(bin).U2)
        conditional_means(bin).U2_mean = mean(bin_data(bin).U2, 2, 'omitnan');
        conditional_means(bin).V2_mean = mean(bin_data(bin).V2, 2, 'omitnan');
    end
end

fprintf('\nComputed conditional mean velocity profiles for all bins.\n');

%% Plotting conditional mean profiles for each bin

figure;
bin = 6; % example bin to plot
plot(conditional_means(bin).U1_mean, Y_profile_air, 'r-', 'LineWidth', 2);
hold on;
plot(conditional_means(bin).U2_mean, Y_profile_liquid, 'b-', 'LineWidth', 2);
ylabel('Y Position (mm)');
xlabel('Mean Velocity Magnitude');
title('conditonal average bin 6');
legend('Air Phase', 'Liquid Phase');
grid on;