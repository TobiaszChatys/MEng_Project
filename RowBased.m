%% import data
clc; clear; close all;

[S, filename] = loadData('L8_G7.mat'); 
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

[~, ~, ~, ~, ~, ~, Y2, ~, ~, ~, ~, ~] = getData(S, Sample_frame);

yrows_liquid = unique(Y2(:, 1));

liquid_count = sum(yrows_liquid >= min_film_height & yrows_liquid <= max_film_height);

fprintf('Number of rows in the liquid region within the film height range: %d\n', liquid_count);
fprintf('Total number of rows within the film height range: %d\n', liquid_count);


% since we have 31 vectors inbetween the min and max film heights, we can determine the number of bins

%% Binning

number_of_vecors_in_bin = 2;
number_of_bins = floor(liquid_count / number_of_vecors_in_bin);
fprintf('Number of bins based on %d vectors per bin: %d\n', number_of_vecors_in_bin, number_of_bins);
bin_edges = linspace(min_film_height, max_film_height, number_of_bins + 1);
fprintf('Bin edges (mm):\n'); disp(bin_edges');
%% plottin

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

%% populate bins with vertical profiles based on film height (parallel)

tic

% Temporary per-frame, per-bin storage for parfor (cells are parfor-safe)
temp_U1 = cell(frames, number_of_bins);
temp_V1 = cell(frames, number_of_bins);
temp_U2 = cell(frames, number_of_bins);
temp_V2 = cell(frames, number_of_bins);

parfor frame = 1:frames
    fprintf('Processing frame %d/%d...\n', frame, frames);
    [X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);

    X_air_columns = X1(1, :);
    X_liquid_columns = X2(1, :);

    % per-frame bin containers
    frame_U1 = cell(1, number_of_bins);
    frame_V1 = cell(1, number_of_bins);
    frame_U2 = cell(1, number_of_bins);
    frame_V2 = cell(1, number_of_bins);

    % air phase
    for column = 1:numel(X_air_columns)
        x_pos = X_air_columns(column);
        local_film_height = interp1(X3, Y3, x_pos, 'linear', 'extrap');
        if isnan(local_film_height)
            continue;
        end
        bin_index = find(local_film_height >= bin_edges(1:end-1) & local_film_height < bin_edges(2:end), 1);
        if isempty(bin_index) && local_film_height == bin_edges(end)
            bin_index = number_of_bins; % assign to last bin if equal to max edge
        elseif isempty(bin_index)
            continue;
        end

        U1_column = U1(:, column);
        V1_column = V1(:, column);

        if isempty(frame_U1{bin_index})
            frame_U1{bin_index} = U1_column;
            frame_V1{bin_index} = V1_column;
        else
            frame_U1{bin_index} = [frame_U1{bin_index}, U1_column];
            frame_V1{bin_index} = [frame_V1{bin_index}, V1_column];
        end
    end

    % liquid phase
    for column = 1:numel(X_liquid_columns)
        x_pos = X_liquid_columns(column);
        local_film_height = interp1(X3, Y3, x_pos, 'linear', 'extrap');
        if isnan(local_film_height)
            continue;
        end
        bin_index = find(local_film_height >= bin_edges(1:end-1) & local_film_height < bin_edges(2:end), 1);
        if isempty(bin_index) && local_film_height == bin_edges(end)
            bin_index = number_of_bins; % assign to last bin if equal to max edge
        elseif isempty(bin_index)
            continue;
        end

        U2_column = U2(:, column);
        V2_column = V2(:, column);

        if isempty(frame_U2{bin_index})
            frame_U2{bin_index} = U2_column;
            frame_V2{bin_index} = V2_column;
        else
            frame_U2{bin_index} = [frame_U2{bin_index}, U2_column];
            frame_V2{bin_index} = [frame_V2{bin_index}, V2_column];
        end
    end

    % write per-frame results into temp arrays (sliced by frame)
    for b = 1:number_of_bins
        temp_U1{frame, b} = frame_U1{b};
        temp_V1{frame, b} = frame_V1{b};
        temp_U2{frame, b} = frame_U2{b};
        temp_V2{frame, b} = frame_V2{b};
    end
end

% aggregate per-frame results into bin_data
for b = 1:number_of_bins

    fprintf('Aggregating data for bin %d/%d...\n', b, number_of_bins);
    for f = 1:frames
        if ~isempty(temp_U1{f, b})
            bin_data(b).U1 = [bin_data(b).U1, temp_U1{f, b}];
        end
        if ~isempty(temp_V1{f, b})
            bin_data(b).V1 = [bin_data(b).V1, temp_V1{f, b}];
        end
        if ~isempty(temp_U2{f, b})
            bin_data(b).U2 = [bin_data(b).U2, temp_U2{f, b}];
        end
        if ~isempty(temp_V2{f, b})
            bin_data(b).V2 = [bin_data(b).V2, temp_V2{f, b}];
        end
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

%% Calculate velocity fluctuations for each bin

fluctuation_data = repmat(struct( ...
    "U1_prime", [], ...
    "V1_prime", [], ...
    "U2_prime", [], ...
    "V2_prime", [] ...
), number_of_bins, 1);

fprintf('\nCalculating velocity fluctuations for each bin...\n');

for bin = 1:number_of_bins
    % Air phase fluctuations
    if ~isempty(bin_data(bin).U1)
        % u' = u_inst - u_mean (subtract mean from each instantaneous profile)
        fluctuation_data(bin).U1_prime = bin_data(bin).U1 - conditional_means(bin).U1_mean;
        fluctuation_data(bin).V1_prime = bin_data(bin).V1 - conditional_means(bin).V1_mean;
        
        fprintf('Bin %d Air: %d fluctuation profiles calculated\n', bin, size(bin_data(bin).U1, 2));
    end
    
    % Liquid phase fluctuations
    if ~isempty(bin_data(bin).U2)
        % u' = u_inst - u_mean (subtract mean from each instantaneous profile)
        fluctuation_data(bin).U2_prime = bin_data(bin).U2 - conditional_means(bin).U2_mean;
        fluctuation_data(bin).V2_prime = bin_data(bin).V2 - conditional_means(bin).V2_mean;
        
        fprintf('Bin %d Liquid: %d fluctuation profiles calculated\n', bin, size(bin_data(bin).U2, 2));
    end
end
%% Calculate RMS (root mean square) of fluctuations

rms_fluctuations = repmat(struct( ...
    "U1_rms", [], ...
    "V1_rms", [], ...
    "U2_rms", [], ...
    "V2_rms", [] ...
), number_of_bins, 1);

for bin = 1:number_of_bins
    % Air phase RMS
    if ~isempty(fluctuation_data(bin).U1_prime)
        rms_fluctuations(bin).U1_rms = sqrt(mean(fluctuation_data(bin).U1_prime.^2, 2, 'omitnan'));
        rms_fluctuations(bin).V1_rms = sqrt(mean(fluctuation_data(bin).V1_prime.^2, 2, 'omitnan'));
    end
    
    % Liquid phase RMS
    if ~isempty(fluctuation_data(bin).U2_prime)
        rms_fluctuations(bin).U2_rms = sqrt(mean(fluctuation_data(bin).U2_prime.^2, 2, 'omitnan'));
        rms_fluctuations(bin).V2_rms = sqrt(mean(fluctuation_data(bin).V2_prime.^2, 2, 'omitnan'));
    end
end

fprintf('\nComputed RMS fluctuations for all bins.\n');

%% Define plotting parameters and labels

color_air = [235, 111, 146] / 255;   % #eb6f92 for air phase
color_liquid = [49, 116, 143] / 255; % #31748f for liquid phase

bin_labels = cell(1, number_of_bins);
for bin = 1:number_of_bins
    bin_labels{bin} = sprintf('Bin %d: %.2f - %.2f mm', bin, bin_edges(bin), bin_edges(bin+1)); 
end

markers = {'o', '+', '*', '.', 'x', 'square', 'diamond', '^', 'p', 'h', 'v', '<', '>', 's', 'd'}; % Define a marker for each bin

%% Plotting mean velocity and RMS profiles for all bins overlaid (log scale)

figure,

subplot(1, 2, 1);
hold on; 

% plot(nan, nan, 'o', 'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', 'none', ...
%     'LineStyle', 'none', 'DisplayName', 'Liquid Phase (hollow)');
% plot(nan, nan, 'o', 'MarkerEdgeColor', color_air, 'MarkerFaceColor', color_air, ...
%     'LineStyle', 'none', 'DisplayName', 'Air Phase (filled)');

% Plot liquid phase for all bins (hollow markers)

for bin = 1:number_of_bins
    if ~isempty(conditional_means(bin).U2_mean)
    semilogx(conditional_means(bin).U2_mean, Y_profile_liquid, 'Marker', ...
        markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_liquid, ...
        'MarkerFaceColor', 'none', 'DisplayName', bin_labels{bin});
    end
end
% Plot air phase for all bins (filled markers)
for bin = 1:number_of_bins
    if ~isempty(conditional_means(bin).U1_mean) 
    semilogx(conditional_means(bin).U1_mean, Y_profile_air, 'Marker', ...
        markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_air, ...
        'MarkerFaceColor', color_air, 'HandleVisibility', 'off');
    end
end

set(gca, 'XScale', 'log');
ylabel('Y Position (mm)', 'FontSize', 12, 'Interpreter', 'latex');
xlabel('$\overline{u}_{(x_0,y)}$ (ms$^{-1}$)', 'FontSize', 12, 'Interpreter', 'latex');
ylim([0, 28]); 
yticks(0:2:28);
xlim([0, 25]);
xticks([0.1 0.5 1 5 10 25]);
title(sprintf('Mean Velocity Profiles With Bins Containing %d Vectors Each', number_of_vecors_in_bin), 'FontSize', 14);
lgd = legend('FontSize', 9); 
lgd.Position(1:2) = [0.32, 0.66];  
grid on;
hold off;

axes('Position', [0.18, 0.5, 0.15, 0.3]); % [left, bottom, width, height]
box on;
hold on;
for bin = 1:number_of_bins
    if ~isempty(conditional_means(bin).U2_mean)
        semilogx(conditional_means(bin).U2_mean, Y_profile_liquid, 'Marker', ...
            markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_liquid, ...
            'MarkerFaceColor', 'none');
    end
end

for bin = 1:number_of_bins
    if ~isempty(conditional_means(bin).U1_mean)
        semilogx(conditional_means(bin).U1_mean, Y_profile_air, 'Marker', ...
            markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_air, ...
            'MarkerFaceColor', color_air);
    end
end
set(gca, 'XScale', 'log');
xlim([0.1, 1]); % Focus on liquid phase velocity range
ylim([min_film_height, max_film_height]); % Focus on film height range
grid on;
xlabel('$\overline{u}$ (ms$^{-1}$)', 'FontSize', 9, 'Interpreter', 'latex');
ylabel('Y (mm)', 'FontSize', 9, 'Interpreter', 'latex');
title('Liquid Phase Detail', 'FontSize', 9);
hold off;

% Plot RMS of fluctuations for liquid phase and air phase for all bins (log scale)
subplot(1, 2, 2);
hold on;

% Plot liquid phase RMS fluctuations (filled markers)
for bin = 1:number_of_bins
    if ~isempty(rms_fluctuations(bin).U2_rms)
    semilogx(rms_fluctuations(bin).U2_rms, Y_profile_liquid, ...
        'Color', color_liquid, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', "none", ...
        'DisplayName', bin_labels{bin});
    end
end

% Plot air phase RMS fluctuations (hollow markers)
for bin = 1:number_of_bins
    if ~isempty(rms_fluctuations(bin).U1_rms) 
    semilogx(rms_fluctuations(bin).U1_rms, Y_profile_air, ...
        'Color', color_air, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'MarkerEdgeColor', color_air, 'MarkerFaceColor', color_air, ...
        'HandleVisibility', 'off');
    end
end

set(gca, 'XScale', 'log');
ylabel('Y Position (mm)', 'FontSize', 12,'Interpreter', 'latex');
xlabel('$u''_{(x_0,y),rms}$ (ms$^{-1}$)', 'FontSize', 12, 'Interpreter', 'latex');
ylim([0, 28]); 
yticks(0:2:28);
xlim([0, 15]);
xticks([0.1 0.5 1 5 10 15]);
title(sprintf('RMS Fluctuations With Bins Containing %d Vectors Each', number_of_vecors_in_bin), 'FontSize', 14);
lgd = legend('FontSize', 9); 
lgd.Position(1:2) = [0.8, 0.7]; 
grid on;
hold off;

% Create inset zoom for RMS subplot
axes('Position', [0.60, 0.5, 0.15, 0.3]); % [left, bottom, width, height]
box on;
hold on;
for bin = 1:number_of_bins
    if ~isempty(rms_fluctuations(bin).U2_rms)
        semilogx(rms_fluctuations(bin).U2_rms, Y_profile_liquid, 'Marker', ...
            markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_liquid, ...
            'MarkerFaceColor', 'none');
    end
end
for bin = 1:number_of_bins
    if ~isempty(rms_fluctuations(bin).U1_rms)
        semilogx(rms_fluctuations(bin).U1_rms, Y_profile_air, 'Marker', ...
            markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_air, ...
            'MarkerFaceColor', color_air);
    end
end
set(gca, 'XScale', 'log');
xlim([0.1, 0.5]); % Focus on liquid phase RMS range
ylim([min_film_height, max_film_height]); % Focus on film height range
grid on;
xlabel('$u''_{rms}$ (ms$^{-1}$)', 'FontSize', 9, 'Interpreter', 'latex');
ylabel('Y (mm)', 'FontSize', 9, 'Interpreter', 'latex');
title('Liquid Phase Detail', 'FontSize', 9);
hold off;

%% Identify representative bins for simplified plotting

min_bin = 1; % Minimum film height (thinnest)
max_bin = number_of_bins; % Maximum film height (thickest)
median_bin = ceil(number_of_bins / 2); % Median film height (middle)

fprintf('\nRepresentative Bins for Simplified Plotting:\n');
fprintf('Min Bin %d: %.2f - %.2f mm\n', min_bin, bin_edges(min_bin), bin_edges(min_bin+1));
fprintf('Median Bin %d: %.2f - %.2f mm\n', median_bin, bin_edges(median_bin), bin_edges(median_bin+1));
fprintf('Max Bin %d: %.2f - %.2f mm\n', max_bin, bin_edges(max_bin), bin_edges(max_bin+1));

%% Plotting for the mix median and max

figure,

subplot(1, 2, 1);
hold on;

% bin_labels = cell(1, number_of_bins);
% for bin = 1:number_of_bins
%     bin_labels{bin} = sprintf('Bin %d: %.2f - %.2f mm', bin, bin_edges(bin), bin_edges(bin+1)); 
% end

Analysis_bins = [min_bin, median_bin, max_bin];
number_of_analysis_bins = length(Analysis_bins);
analysis_bin_labels = cell(1, number_of_analysis_bins);

for index = 1:number_of_analysis_bins
    bin = Analysis_bins(index);
    if bin == min_bin
        label_type = 'min';
    elseif bin == median_bin
        label_type = 'median';
    else
        label_type = 'max';
    end
    analysis_bin_labels{index} = sprintf('Bin %d: %.2f - %.2f mm (%s)', bin, bin_edges(bin), bin_edges(bin+1), label_type);
end

% Plot liquid phase for all bins (hollow markers)

for index = 1:number_of_analysis_bins
    bin = Analysis_bins(index);
    if ~isempty(conditional_means(bin).U2_mean)
    semilogx(conditional_means(bin).U2_mean, Y_profile_liquid, 'Marker', ...
        markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_liquid, ...
        'MarkerFaceColor', 'none', 'DisplayName', analysis_bin_labels{index});
    end
end

% Plot air phase for all bins (filled markers)

for index = 1:number_of_analysis_bins
    bin = Analysis_bins(index);
    if ~isempty(conditional_means(bin).U1_mean)
    semilogx(conditional_means(bin).U1_mean, Y_profile_air, 'Marker', ...
        markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_air, ...
        'MarkerFaceColor', color_air, 'HandleVisibility', 'off');
    end
end

set(gca, 'XScale', 'log');
ylabel('Y Position (mm)', 'FontSize', 12, 'Interpreter', 'latex');
xlabel('$\overline{u}_{(x_0,y)}$ (ms$^{-1}$)', 'FontSize', 12, 'Interpreter', 'latex');
ylim([0, 28]); 
yticks(0:2:28);
xlim([0, 25]);
xticks([0.1 0.5 1 5 10 25]);
title(sprintf('Mean Velocity Profiles With Bins Containing %d Vectors Each', number_of_vecors_in_bin), 'FontSize', 14);
lgd = legend('FontSize', 9); 
lgd.Position(1:2) = [0.32, 0.66];  
grid on;
hold off;

axes('Position', [0.18, 0.5, 0.15, 0.3]); % [left, bottom, width, height]
box on;
hold on;

for index = 1:number_of_analysis_bins
    bin = Analysis_bins(index);
    if ~isempty(conditional_means(bin).U2_mean)
        semilogx(conditional_means(bin).U2_mean, Y_profile_liquid, 'Marker', ...
            markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_liquid, ...
            'MarkerFaceColor', 'none');
    end
end

for index = 1:number_of_analysis_bins
    bin = Analysis_bins(index);
    if ~isempty(conditional_means(bin).U1_mean)
        semilogx(conditional_means(bin).U1_mean, Y_profile_air, 'Marker', ...
            markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_air, ...
            'MarkerFaceColor', color_air);
    end
end

set(gca, 'XScale', 'log');
xlim([0.1, 1]); % Focus on liquid phase velocity range
ylim([min_film_height, max_film_height]); % Focus on film height range
grid on;
xlabel('$\overline{u}$ (ms$^{-1}$)', 'FontSize', 9, 'Interpreter', 'latex');
ylabel('Y (mm)', 'FontSize', 9, 'Interpreter', 'latex');
title('Liquid Phase Detail', 'FontSize', 9);
hold off;

subplot(1, 2, 2);
hold on;

% Plot liquid phase RMS fluctuations (filled markers)
for index = 1:number_of_analysis_bins
bin = Analysis_bins(index);
    if ~isempty(rms_fluctuations(bin).U2_rms)
    semilogx(rms_fluctuations(bin).U2_rms, Y_profile_liquid, ...
        'Color', color_liquid, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', "none", ...
        'DisplayName', analysis_bin_labels{index});
    end
end

% Plot air phase RMS fluctuations (hollow markers)
for index = 1:number_of_analysis_bins
    bin = Analysis_bins(index);
    if ~isempty(rms_fluctuations(bin).U1_rms) 
    semilogx(rms_fluctuations(bin).U1_rms, Y_profile_air, ...
        'Color', color_air, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'MarkerEdgeColor', color_air, 'MarkerFaceColor', color_air, ...
        'HandleVisibility', 'off');
    end
end

set(gca, 'XScale', 'log');
ylabel('Y Position (mm)', 'FontSize', 12,'Interpreter', 'latex');
xlabel('$u''_{(x_0,y),rms}$ (ms$^{-1}$)', 'FontSize', 12, 'Interpreter', 'latex');
ylim([0, 28]); 
yticks(0:2:28);
xlim([0, 10]);
xticks([0.1 0.5 1 5 10 50 100]);
title(sprintf('RMS Fluctuations With Bins Containing %d Vectors Each', number_of_vecors_in_bin), 'FontSize', 14);
lgd = legend('FontSize', 9); 
lgd.Position(1:2) = [0.8, 0.7];
grid on;
hold off;

% Create inset zoom for RMS subplot
axes('Position', [0.60, 0.5, 0.15, 0.3]); % [left, bottom, width, height]
box on;
hold on;
for index = 1:number_of_analysis_bins
    bin = Analysis_bins(index);
    if ~isempty(rms_fluctuations(bin).U2_rms)
        semilogx(rms_fluctuations(bin).U2_rms, Y_profile_liquid, 'Marker', ...
            markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_liquid, ...
            'MarkerFaceColor', 'none');
    end
end
for index = 1:number_of_analysis_bins
    bin = Analysis_bins(index);
    if ~isempty(rms_fluctuations(bin).U1_rms)
        semilogx(rms_fluctuations(bin).U1_rms, Y_profile_air, 'Marker', ...
            markers{bin}, 'LineStyle', 'none', 'MarkerEdgeColor', color_air, ...
            'MarkerFaceColor', color_air);
    end
end
set(gca, 'XScale', 'log');
xlim([0.1, 0.5]); % Focus on liquid phase RMS range
ylim([min_film_height, max_film_height]); % Focus on film height range
grid on;
xlabel('$u''_{rms}$ (ms$^{-1}$)', 'FontSize', 9, 'Interpreter', 'latex');
ylabel('Y (mm)', 'FontSize', 9, 'Interpreter', 'latex');
title('Liquid Phase Detail', 'FontSize', 9);
hold off;
