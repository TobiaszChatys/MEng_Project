%% import data
clc; clear; close all;

[S, filename] = loadData('L8_G3.mat'); 
frames = size(S.all_u_matrix_liquid, 3);


%% compute Film height stats

Film_height_matrix = S.smoothed_film_height_matrix_out * 1e3;
all_film_heights = Film_height_matrix(Film_height_matrix >= 0); % Filter out negative heights
invalid_film_heights = Film_height_matrix(Film_height_matrix < 0); % Count negative heights as invalid


min_film_height = min(all_film_heights);
max_film_height = max(all_film_heights);
mean_film_height = mean(all_film_heights);

fprintf('Film Height Analysis for %s:\n', filename);
fprintf('We are expecting 4932500 film height measurements (1973 x 2500) across all frames.\n');
fprintf('Total Number of Film Height Measurements: %d\n', numel(all_film_heights));
fprintf('Number of Invalid Film Height Measurements: %d\n', numel(invalid_film_heights));
fprintf('Total number of film heights after filtering out invalid measurements: %d\n', numel(all_film_heights) + numel(invalid_film_heights));
fprintf('Film Height Statistics:\n');
fprintf('Minimum Film Height: %.2f mm\n', min_film_height);
fprintf('Maximum Film Height: %.2f mm\n', max_film_height);
fprintf('Mean Film Height: %.2f mm\n', mean_film_height);

%% Velocity Vectors in range of film heights

Sample_frame = 40;

[~, Y1, ~, ~, ~, ~, Y2, ~, ~, ~, ~, Y3] = getData(S, Sample_frame);

yrows_liquid = unique(Y2(:, 1));

liquid_count = sum(yrows_liquid >= min_film_height & yrows_liquid <= max_film_height);

fprintf('Number of rows in the liquid region within the film height range: %d\n', liquid_count);
fprintf('Total number of rows within the film height range: %d\n', liquid_count);


% since we have 31 vectors inbetween the min and max film heights, we can determine the number of bins

%% Binning
bin_edges = linspace(min_film_height, max_film_height, 8);
fprintf('Bin edges for film height range: %s\n', mat2str(bin_edges));

number_of_bins = length(bin_edges) - 1;

for bin = 1:number_of_bins
    bin_counts(bin) = sum(all_film_heights >= bin_edges(bin) & all_film_heights < bin_edges(bin+1));
    fprintf('Bin %d (%.2f - %.2f mm): %d data points\n', bin, bin_edges(bin), bin_edges(bin+1), bin_counts(bin));
end



%% plotting

frame = 10;
[X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);

quiver(X1 + 11, Y1, U1, V1, 0.8, 'k');
hold on;
quiver(X2 + 11, Y2, U2, V2, 0.8, 'b');
hold on; 
plot(X3 + 11, Y3, 'g', 'LineWidth', 4)
hold off;
yline(min_film_height, 'r--', 'LineWidth', 2);
yline(max_film_height, 'r--', 'LineWidth', 2);

for bin = 2:length(bin_edges)-1
    yline(bin_edges(bin), 'm--', 'LineWidth', 1);
end

xlabel('X Position (mm)');
ylabel('Y Position (mm)');
title(sprintf('L8-G3 Velocity Vector Map - Frame %d', frame));
xlim([-15 15]);
ylim([0 28]);

y_ticks = 0:2:28;
yticks(y_ticks);
yticklabels(y_ticks)
xline(0, 'k--', 'LineWidth', 2);

%% Storage for binning results

temp_U1_cell = cell(1, number_of_bins);
temp_V1_cell = cell(1, number_of_bins);
temp_U2_cell = cell(1, number_of_bins);
temp_V2_cell = cell(1, number_of_bins);

bin_edges_local = bin_edges; % Store bin edges in a local variable for use in parfor loop
number_of_bins_local = number_of_bins; % Store number of bins in a local variable for use in parfor loop
%% Binning velocity vectors based on film height

tic
parfor frame = 1:frames

    [X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);
    
    x_air_columns = X1(1, :);
    y_liquid_columns = Y2(1, :);

    % Temporary storage for current frame's binned data
    frame_U1_bins = cell(number_of_bins_local, 1);
    frame_V1_bins = cell(number_of_bins_local, 1);


    % Air phase 

    for column = 1:numel(x_air_columns)
        x_position = x_air_columns(column);
        film_height_at_column = interp1(X3, Y3, x_position, 'linear', 'extrap'); % Get film height at this x position
        bin_index = find(film_height_at_column >= bin_edges_local(1:end-1) & film_height_at_column < bin_edges_local(2:end), 1);
        if ~isempty(bin_index)
            frame_U1_bins{bin_index} = [frame_U1_bins{bin_index}; U1(:, column)];
            frame_V1_bins{bin_index} = [frame_V1_bins{bin_index}; V1(:, column)];
        end
    end

    for b = 1:number_of_bins_local
        temp_U1{frame, b} = frame_U1_bins{b};
        temp_V1{frame, b} = frame_V1_bins{b};
    end 

end
toc
