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
h_std = std(heights);

fprintf('Film Height Statistics (All %d Frames):\n', frames);
fprintf('Total Data Points: %d\n', length(heights));
fprintf('Max Height: %.2f mm\n', h_max);
fprintf('Min Height: %.2f mm\n', h_min);
fprintf('Mean Height: %.2f mm\n', h_mean);
fprintf('Std Dev of Height: %.2f mm\n', h_std);

%% Create bins for conditional averaging

clc

bin_edges = [h_min, 1.0, 1.5, 2.0, 2.5, 3.5, h_max];
n_bins = length(bin_edges) - 1;

% count data points in each bin
bin_counts = zeros(n_bins, 1);
fprintf('\nData Points in Each Bin:\n');

for bin = 1:n_bins
    bin_counts(bin) = sum(heights >= bin_edges(bin) & heights < bin_edges(bin+1));
    fprintf('Bin %d (%.2f - %.2f mm): %d data points\n', bin, bin_edges(bin), bin_edges(bin+1), bin_counts(bin));
end


%% Visualising height distribution
figure;
histogram(heights, 50);
xlabel('Film Height (mm)');
ylabel('Frequency');
title(sprintf('Distribution of Film Heights (All %d Frames)', frames));
grid on;

% Add lines for mean, max, min
hold on;
xline(h_mean, 'r--', 'Mean', 'LineWidth', 1);
xline(h_max, 'g--', 'Max', 'LineWidth', 1);
xline(h_min, 'b--', 'Min', 'LineWidth', 1);
hold off;

% add bin edges
for bin = 2:length(bin_edges)-1
    xline(bin_edges(bin), 'k:', 'LineWidth', 1.5);
end

%% Create Storage for velocity fields in each bin

bin_data = struct();

for bin = 1:n_bins
    bin_data(bin).U1 = [];
    bin_data(bin).V1 = [];
    bin_data(bin).X1 = [];
    bin_data(bin).Y1 = [];
    
    bin_data(bin).U2 = [];
    bin_data(bin).V2 = [];
    bin_data(bin).X2 = [];
    bin_data(bin).Y2 = [];

end

%% Populate bins with velocity data based on height conditions

for frame = 1:10
    [U1, V1, X1, Y1, U2, V2, X2, Y2, ~, ~, ~, ~] = getData(S, frame);

    % for each spatial point in the liquid phase 
    for i = 1:length(X1(:))
        % get the height at this point
        h = Y1(i);
        for bin = 1:n_bins
            if h >= bin_edges(bin) && h < bin_edges(bin+1)
                bin_data(bin).U1 = [bin_data(bin).U1; U1(i)];
                bin_data(bin).V1 = [bin_data(bin).V1; V1(i)];
                bin_data(bin).X1 = [bin_data(bin).X1; X1(i)];
                bin_data(bin).Y1 = [bin_data(bin).Y1; Y1(i)];
                break;
            end
        end
    end

    % for each spatial point in the gas phase
    for i = 1:length(X2(:))
        % get the height at this point
        h = Y2(i);
        for bin = 1:n_bins
            if h >= bin_edges(bin) && h < bin_edges(bin+1)
                bin_data(bin).U2 = [bin_data(bin).U2; U2(i)];
                bin_data(bin).V2 = [bin_data(bin).V2; V2(i)];
                bin_data(bin).X2 = [bin_data(bin).X2; X2(i)];
                bin_data(bin).Y2 = [bin_data(bin).Y2; Y2(i)];
                break;
            end
        end
    end
end

%% Verify data collection in bins

fprintf('\nVelocity Data Points Collected in Each Bin:\n');
for bin = 1:n_bins
    fprintf('Bin %d (%.2f - %.2f mm):\n', bin, bin_edges(bin), bin_edges(bin+1));
    fprintf('  Liquid points: %d\n', length(bin_data(bin).X1));
    fprintf('  Gas points: %d\n', length(bin_data(bin).X2));
end