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

xlabel('X Position (mm)');
ylabel('Y Position (mm)');
title(sprintf('L8-G3 Velocity Vector Map - Frame %d', frame));
xlim([-15 15]);
ylim([0 28]);

y_ticks = 0:2:28;
yticks(y_ticks);
yticklabels(y_ticks)
xline(0, 'k--', 'LineWidth', 2);

%% Velocity Vectors in range of film heights

Sample_frame = 10;

[~, Y1, ~, ~, ~, ~, Y2, ~, ~, ~, ~, Y3] = getData(S, Sample_frame);

yrows_air = unique(Y1(:, 1));
yrows_liquid = unique(Y2(:, 1));

air_count = sum(yrows_air >= min_film_height & yrows_air <= max_film_height);
liquid_count = sum(yrows_liquid >= min_film_height & yrows_liquid <= max_film_height);
total_rows = air_count + liquid_count;
fprintf('Number of rows in the air region within the film height range: %d\n', air_count);
fprintf('Number of rows in the liquid region within the film height range: %d\n', liquid_count);
fprintf('Total number of rows within the film height range: %d\n', total_rows);


