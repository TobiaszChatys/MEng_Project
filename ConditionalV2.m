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

% Define bins using percentiles
film_p25 = prctile(all_film_heights, 25);
film_p75 = prctile(all_film_heights, 75);
film_median = median(all_film_heights);

% Create bin edges using percentiles as outer boundaries
bin_edges = [
    film_min_height,  % minimum
    film_p25,         % 25th percentile
    film_median,      % 50th percentile (median)
    film_p75,         % 75th percentile
    film_max_height   % maximum
];

number_of_bins = length(bin_edges) - 1;

fprintf('\nPercentile-based Bin Edges:\n');
fprintf('25th percentile: %.2f mm\n', film_p25);
fprintf('50th percentile (median): %.2f mm\n', film_median);
fprintf('75th percentile: %.2f mm\n', film_p75);
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
bin = 1; % example bin to plot
plot(conditional_means(bin).U1_mean, Y_profile_air, 'r-', 'LineWidth', 2);
hold on;
plot(conditional_means(bin).U2_mean, Y_profile_liquid, 'b-', 'LineWidth', 2);
ylabel('Y Position (mm)');
xlabel('Mean Velocity Magnitude');
title('conditonal average bin 1');
legend('Air Phase', 'Liquid Phase');
grid on;

%% Define plotting styles for each bin

% Same markers for both phases
markers = {'s', 'o', '^', 'd'};

% Rose Pine colors
% Each bin gets one color used for both air (hollow) and liquid (filled)
% colors = [
%     0.85, 0.52, 0.60;  % rose (darker)
%     0.90, 0.77, 0.53;  % gold (darker)
%     0.76, 0.56, 0.78;  % iris
%     0.56, 0.73, 0.78   % foam
% ];

% Black for air, Blue for liquid
color_air = [215, 130, 126] / 255;   % #d7827e for air phase
color_liquid = [40, 105, 131] / 255;   % #286983 for liquid phase
color_single = [180,  99, 122] / 255;  % #b4637a for single-height profile


percentile_labels = {
    sprintf('%.2f - %.2f mm', bin_edges(1), bin_edges(2)), ...
    sprintf('%.2f - %.2f mm', bin_edges(2), bin_edges(3)), ...
    sprintf('%.2f - %.2f mm', bin_edges(3), bin_edges(4)), ...
    sprintf('%.2f - %.2f mm', bin_edges(4), bin_edges(5))
};

%% Plotting all bins overlaid - LOG SCALE

figure('Position', [100, 100, 900, 600]);
hold on;

% Plot liquid phase for all bins (hollow markers)

for bin = 1:number_of_bins
    if ~isempty(conditional_means(bin).U2_mean)
    semilogx(conditional_means(bin).U2_mean, Y_profile_liquid, ...
        'Color', color_liquid, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 1, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', 'none', ...
        'DisplayName', percentile_labels{bin});
    end
end

% Plot air phase for all bins (filled markers)
for bin = 1:number_of_bins
    if ~isempty(conditional_means(bin).U1_mean) 
    semilogx(conditional_means(bin).U1_mean, Y_profile_air, ...
        'Color', color_air, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_air, 'MarkerFaceColor', color_air, ...
        'HandleVisibility', 'off');
    end
end


% horizontal reference line at bin edges

for bin = 2:number_of_bins
    yline(bin_edges(bin), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
set(gca, 'XScale', 'log');
ylabel('Y Position (mm)', 'FontSize', 12);
xlabel('Mean Velocity Magnitude (log scale)', 'FontSize', 12);
ylim([0, 28]); 
yticks(0:2:28);
xticks([0.1 0.5 1 5 10 50 100]);
title('Conditional Averaged Velocity Profiles (Log Scale) - All Bins Overlaid', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

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

%% Plot RMS velocity fluctuations - LOG SCALE

figure('Position', [100, 100, 900, 600]);
hold on;

% Plot air phase RMS fluctuations (hollow markers)
for bin = 1:number_of_bins
    if ~isempty(rms_fluctuations(bin).U1_rms) 
    semilogx(rms_fluctuations(bin).U1_rms, Y_profile_air, ...
        'Color', color_air, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_air, 'MarkerFaceColor', 'none', ...
        'HandleVisibility', 'off');
    end
end

% Plot liquid phase RMS fluctuations (filled markers)
for bin = 1:number_of_bins
    if ~isempty(rms_fluctuations(bin).U2_rms)
    semilogx(rms_fluctuations(bin).U2_rms, Y_profile_liquid, ...
        'Color', color_liquid, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 1, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', color_liquid, ...
        'DisplayName', percentile_labels{bin});
    end
end

set(gca, 'XScale', 'log');
ylabel('Y Position (mm)', 'FontSize', 12);
xlabel('RMS Velocity Fluctuation (log scale)', 'FontSize', 12);
ylim([0, 28]); 
yticks(0:2:28);
xticks([0.1 0.5 1 5 10 50 100]);
title('RMS Velocity Fluctuations - All Bins Overlaid', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

%% Calculate turbulence intensity for each bin

turbulence_intensity = repmat(struct( ...
    "U1_TI", [], ...
    "V1_TI", [], ...
    "U2_TI", [], ...
    "V2_TI", [] ...
), number_of_bins, 1);

fprintf('\nCalculating turbulence intensity for each bin...\n');

for bin = 1:number_of_bins
    % Air phase turbulence intensity
    if ~isempty(rms_fluctuations(bin).U1_rms)
        % TI = RMS / mean velocity
        turbulence_intensity(bin).U1_TI = rms_fluctuations(bin).U1_rms ./ abs(conditional_means(bin).U1_mean);
        turbulence_intensity(bin).V1_TI = rms_fluctuations(bin).V1_rms ./ abs(conditional_means(bin).V1_mean);
        
        fprintf('Bin %d Air: Turbulence intensity calculated\n', bin);
    end
    
    % Liquid phase turbulence intensity
    if ~isempty(rms_fluctuations(bin).U2_rms)
        % TI = RMS / mean velocity
        turbulence_intensity(bin).U2_TI = rms_fluctuations(bin).U2_rms ./ abs(conditional_means(bin).U2_mean);
        turbulence_intensity(bin).V2_TI = rms_fluctuations(bin).V2_rms ./ abs(conditional_means(bin).V2_mean);
        
        fprintf('Bin %d Liquid: Turbulence intensity calculated\n', bin);
    end
end

fprintf('\nComputed turbulence intensity for all bins.\n');

%% Calculate single height profile comparison 

target_film_height = 3.0; % mm
fprintf('\nCalculating single-height profile at film height: %.2f mm\n', target_film_height);

% Find all data points within a narrow range around 3mm (e.g., ±0.1 mm)
height_tolerance = 0.1; % mm

single_height_data = struct( ...
    'U1', [], ...
    'V1', [], ...
    'U2', [], ...
    'V2', [] ...
);

% Extract profiles at this specific height across all frames
for frame = 1:frames
    [X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);
    
    X_air_columns = X1(1, :);
    X_liquid_columns = X2(1, :);
    
    % Air phase at 3mm height
    for column = 1:size(U1, 2)
        x_pos = X_air_columns(column);
        local_film_height = interp1(X3, Y3, x_pos, 'linear', 'extrap');
        
        if abs(local_film_height - target_film_height) < height_tolerance
            U1_column = U1(:, column);
            V1_column = V1(:, column);
            single_height_data.U1 = [single_height_data.U1, U1_column];
            single_height_data.V1 = [single_height_data.V1, V1_column];
        end
    end
    
    % Liquid phase at 3mm height
    for column = 1:size(U2, 2)
        x_pos = X_liquid_columns(column);
        local_film_height = interp1(X3, Y3, x_pos, 'linear', 'extrap');
        
        if abs(local_film_height - target_film_height) < height_tolerance
            U2_column = U2(:, column);
            V2_column = V2(:, column);
            single_height_data.U2 = [single_height_data.U2, U2_column];
            single_height_data.V2 = [single_height_data.V2, V2_column];
        end
    end
end

% Compute mean profiles for single height
single_height_means = struct( ...
    'U1_mean', mean(single_height_data.U1, 2, 'omitnan'), ...
    'V1_mean', mean(single_height_data.V1, 2, 'omitnan'), ...
    'U2_mean', mean(single_height_data.U2, 2, 'omitnan'), ...
    'V2_mean', mean(single_height_data.V2, 2, 'omitnan') ...
);

% Compute RMS for single height
single_height_fluctuations = struct( ...
    'U1_prime', single_height_data.U1 - single_height_means.U1_mean, ...
    'V1_prime', single_height_data.V1 - single_height_means.V1_mean, ...
    'U2_prime', single_height_data.U2 - single_height_means.U2_mean, ...
    'V2_prime', single_height_data.V2 - single_height_means.V2_mean ...
);

single_height_rms = struct( ...
    'U1_rms', sqrt(mean(single_height_fluctuations.U1_prime.^2, 2, 'omitnan')), ...
    'V1_rms', sqrt(mean(single_height_fluctuations.V1_prime.^2, 2, 'omitnan')), ...
    'U2_rms', sqrt(mean(single_height_fluctuations.U2_prime.^2, 2, 'omitnan')), ...
    'V2_rms', sqrt(mean(single_height_fluctuations.V2_prime.^2, 2, 'omitnan')) ...
);

single_height_ti = struct( ...
    'U1_ti', single_height_rms.U1_rms ./ abs(single_height_means.U1_mean), ...
    'V1_ti', single_height_rms.V1_rms ./ abs(single_height_means.V1_mean), ...
    'U2_ti', single_height_rms.U2_rms ./ abs(single_height_means.U2_mean), ...
    'V2_ti', single_height_rms.V2_rms ./ abs(single_height_means.V2_mean) ...
);

fprintf('Single-height profile calculated from %d profiles\n', size(single_height_data.U1, 2) + size(single_height_data.U2, 2));
%% Plot turbulence intensity

figure('Position', [100, 100, 900, 600]);
hold on;

% Plot air phase turbulence intensity (hollow markers)
for bin = 1:number_of_bins
    if ~isempty(turbulence_intensity(bin).U1_TI) 
    semilogx(turbulence_intensity(bin).U1_TI, Y_profile_air, ...
        'Color', color_air, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_air, 'MarkerFaceColor', 'none', ...
        'HandleVisibility', 'off');
    end
end

% Plot liquid phase turbulence intensity (filled markers)
for bin = 1:number_of_bins
    if ~isempty(turbulence_intensity(bin).U2_TI)
    semilogx(turbulence_intensity(bin).U2_TI, Y_profile_liquid, ...
        'Color', color_liquid, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 1, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', color_liquid, ...
        'DisplayName', percentile_labels{bin});
    end
end

set(gca, 'XScale', 'log');
ylabel('Y Position (mm)', 'FontSize', 12);
xlabel('Turbulence Intensity (u_{rms} / U_{mean})', 'FontSize', 12);
ylim([0, 28]); 
yticks(0:2:28);
xticks([0.1 0.5 1 5 10 50 100]);
title('Turbulence Intensity - All Bins Overlaid', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;
%% Plot Mean Velocity, RMS Fluctuations, and Turbulence Intensity side by side

figure;

% First subplot: Mean Velocity
subplot(1, 3, 1);
hold on;

% Plot liquid phase for all bins (filled markers)
for bin = 1:number_of_bins
    if ~isempty(conditional_means(bin).U2_mean)
    semilogx(conditional_means(bin).U2_mean, Y_profile_liquid, ...
        'Color', color_liquid, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 1, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', "none", ...
        'DisplayName', percentile_labels{bin});
    end
end

% Plot air phase for all bins (hollow markers)
for bin = 1:number_of_bins
    if ~isempty(conditional_means(bin).U1_mean) 
    semilogx(conditional_means(bin).U1_mean, Y_profile_air, ...
        'Color', color_air, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_air, 'MarkerFaceColor', color_air, ...
        'HandleVisibility', 'off');
    end
end

if ~isempty(single_height_means.U1_mean)
    semilogx(single_height_means.U1_mean, Y_profile_air, ...
        'Color', color_single, 'LineStyle', '-', 'LineWidth', 3, ...
        'DisplayName', sprintf('Single-Height (%.2f mm)', target_film_height));
end

if ~isempty(single_height_means.U2_mean)
    semilogx(single_height_means.U2_mean, Y_profile_liquid, ...
        'Color', color_single, 'LineStyle', '-', 'LineWidth', 3, ...
        'HandleVisibility', 'off');
end

% Add horizontal reference lines
for bin = 2:number_of_bins
    yline(bin_edges(bin), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

set(gca, 'XScale', 'log');
xlim([0.01, 50]);
ylabel('Y Position (mm)', 'FontSize', 12);
xlabel('$\overline{u}_{(x_0,y)}$ (ms$^{-1}$)', 'FontSize', 12, 'Interpreter', 'latex');
ylim([0, 28]); 
yticks(0:2:28);
xticks([0.1 0.5 1 5 10 50]);
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% Second subplot: RMS Fluctuations
subplot(1, 3, 2);
hold on;

% Plot liquid phase RMS fluctuations (filled markers)
for bin = 1:number_of_bins
    if ~isempty(rms_fluctuations(bin).U2_rms)
    semilogx(rms_fluctuations(bin).U2_rms, Y_profile_liquid, ...
        'Color', color_liquid, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 1, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', "none", ...
        'DisplayName', percentile_labels{bin});
    end
end

% Plot air phase RMS fluctuations (hollow markers)
for bin = 1:number_of_bins
    if ~isempty(rms_fluctuations(bin).U1_rms) 
    semilogx(rms_fluctuations(bin).U1_rms, Y_profile_air, ...
        'Color', color_air, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_air, 'MarkerFaceColor', color_air, ...
        'HandleVisibility', 'off');
    end
end

% Overlay single-height RMS as thick lines
if ~isempty(single_height_rms.U1_rms)
    semilogx(single_height_rms.U1_rms, Y_profile_air, ...
        'Color', color_single, 'LineStyle', '-', 'LineWidth', 3, ...
        'DisplayName',sprintf('Single-Height (%.2f mm)', target_film_height));
end

if ~isempty(single_height_rms.U2_rms)
    semilogx(single_height_rms.U2_rms, Y_profile_liquid, ...
        'Color', color_single, 'LineStyle', '-', 'LineWidth', 3, ...
        'HandleVisibility', 'off');
end

% Add horizontal reference lines
for bin = 2:number_of_bins
    yline(bin_edges(bin), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

set(gca, 'XScale', 'log');
xlim([0.01, 10]);
ylabel('Y Position (mm)', 'FontSize', 12);
xlabel('$u''_{(x_0,y),rms}$ (ms$^{-1}$)', 'FontSize', 12, 'Interpreter', 'latex');
ylim([0, 28]); 
yticks(0:2:28);
xticks([0.1 0.5 1 5 10 50]);
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% Third subplot: Turbulence Intensity
subplot(1, 3, 3);
hold on;

% Plot liquid phase turbulence intensity (filled markers)
for bin = 1:number_of_bins
    if ~isempty(turbulence_intensity(bin).U2_TI)
    semilogx(turbulence_intensity(bin).U2_TI, Y_profile_liquid, ...
        'Color', color_liquid, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 1, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', "none", ...
        'DisplayName', percentile_labels{bin});
    end
end

% Plot air phase turbulence intensity (hollow markers)
for bin = 1:number_of_bins
    if ~isempty(turbulence_intensity(bin).U1_TI) 
    semilogx(turbulence_intensity(bin).U1_TI, Y_profile_air, ...
        'Color', color_air, 'Marker', markers{bin}, 'LineStyle', 'none', ...
        'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerEdgeColor', color_air, 'MarkerFaceColor', color_air, ...
        'HandleVisibility', 'off');
    end
end

% Overlay single-height TI as thick lines
if ~isempty(single_height_ti.U1_ti)
    semilogx(single_height_ti.U1_ti, Y_profile_air, ...
        'Color', color_single, 'LineStyle', '-', 'LineWidth', 3, ...
        'DisplayName', sprintf('Single-Height (%.2f mm)', target_film_height));
end

if ~isempty(single_height_ti.U2_ti)
    semilogx(single_height_ti.U2_ti, Y_profile_liquid, ...
        'Color', color_single, 'LineStyle', '-', 'LineWidth', 3, ...
        'HandleVisibility', 'off');
end

% Add horizontal reference lines
for bin = 2:number_of_bins
    yline(bin_edges(bin), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

set(gca, 'XScale', 'log');
xlim([0.01, 10]);
ylabel('Y Position (mm)', 'FontSize', 12);
xlabel('$I = \frac{u''_{(x_0,y),rms}}{\overline{u}_{(x_0,y)}}$', 'FontSize', 12, 'Interpreter', 'latex');
ylim([0, 28]); 
yticks(0:2:28);
xticks([0.1 0.5 1 5 10 50]);
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

