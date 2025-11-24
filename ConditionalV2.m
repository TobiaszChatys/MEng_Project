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

bin_data = struct( ...
    "U1", [], ...
    "V1", [], ...    
    "U2", [], ...
    "V2", [] ...
);

