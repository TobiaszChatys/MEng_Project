clc; clear; close all;

[S, filename] = loadData('L8_G9.mat'); 
frames = size(S.smoothed_film_height_matrix_out, 2);

heights = [];

for frame = 1:frames
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, Y3] = getData(S, frame);
    heights = [heights; Y3(:)];
end


% stats
h_max = max(heights);
h_min = min(heights);
h_mean = mean(heights);
h_std = std(heights);

fprintf('Film Height Statistics at Frame %d:\n', frame);
fprintf('Max Height: %.2f mm\n', h_max);
fprintf('Min Height: %.2f mm\n', h_min);
fprintf('Mean Height: %.2f mm\n', h_mean);
fprintf('Std Dev of Height: %.2f mm\n', h_std);

% Visualising height distribution
figure;
histogram(heights, 5);
xlabel('Film Height (mm)');
ylabel('Frequency');
title('Distribution of Film Heights at Frame 339');
grid on;

% Add lines for mean, max, min
hold on;
xline(h_mean, 'r-', 'Mean', 'LineWidth', 2);
xline(h_max, 'g-', 'Max', 'LineWidth', 2);
xline(h_min, 'b-', 'Min', 'LineWidth', 2);
hold off;

% % Plot Velocity Vectors
% figure
% quiver(X1 + 11, Y1, U1, V1, 0.8, 'k');
% hold on;
% quiver(X2 + 11, Y2, U2, V2, 0.8, 'b');
% hold on;
% plot(X3 + 11, Y3, 'g', 'LineWidth', 4)
% hold off;
% xlim([-15 15]);
% xlabel('X Position (mm)');
% ylabel('Y Position (mm)');
% title('L8-G3 Velocity Vector Map');
% ylim([0 28]);
% y_ticks = 0:2:28;
% yticks(y_ticks);
% yticklabels(y_ticks)
% xline(0, 'k--', 'LineWidth', 2);