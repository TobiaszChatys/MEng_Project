clc; clear; close all;

[S, filename] = loadData('L8_G9.mat');

frames = size(S.all_u_matrix_liquid, 3);


% picking a fixed column to perform simple time-averaging on
column_index = 25;

% limit data to valid index:
U1_air_line = squeeze(S.all_u_matrix_air(:, column_index, 1:frames));
V1_air_line = squeeze(S.all_v_matrix_air(:, column_index, 1:frames));

% Set zero velocities to NaN
zero_mask = (U1_air_line == 0) & (V1_air_line == 0);
U1_air_line(zero_mask) = NaN;
V1_air_line(zero_mask) = NaN;

U2_liquid_line = squeeze(S.all_u_matrix_liquid(:, column_index, 1:frames));
V2_liquid_line = squeeze(S.all_v_matrix_liquid(:, column_index, 1:frames));

% Set zero velocities to NaN
zero_mask = (U2_liquid_line == 0) & (V2_liquid_line == 0);
U2_liquid_line(zero_mask) = NaN;
V2_liquid_line(zero_mask) = NaN;


u1_air_mean = mean(U1_air_line, 2, 'omitnan');
v1_air_mean = mean(V1_air_line, 2, 'omitnan');
u2_liquid_mean = mean(U2_liquid_line, 2, 'omitnan');
v2_liquid_mean = mean(V2_liquid_line, 2, 'omitnan');


mean_air_velocity_magnitude = hypot(u1_air_mean, v1_air_mean);
mean_liquid_velocity_magnitude = hypot(u2_liquid_mean, v2_liquid_mean);

% get y positions
y_air_line = squeeze(S.all_transposed_y_position_matrix_air(:, column_index, 1)) * 1000;
y_liquid_line = squeeze(S.all_transposed_y_position_matrix_liquid(:, column_index, 1)) * 1000;

%% Calculate velocity fluctuations

U1_air_fluctuations = U1_air_line - u1_air_mean;
V1_air_fluctuations = V1_air_line - v1_air_mean;
U2_liquid_fluctuations = U2_liquid_line - u2_liquid_mean;
V2_liquid_fluctuations = V2_liquid_line - v2_liquid_mean;

% Calculate RMS fluctuations
U1_air_rms = sqrt(mean(U1_air_fluctuations.^2, 2, 'omitnan'));
V1_air_rms = sqrt(mean(V1_air_fluctuations.^2, 2, 'omitnan'));
U2_liquid_rms = sqrt(mean(U2_liquid_fluctuations.^2, 2, 'omitnan'));
V2_liquid_rms = sqrt(mean(V2_liquid_fluctuations.^2, 2, 'omitnan'));

% Calculate turbulence intensity
TI_air = U1_air_rms ./ abs(u1_air_mean);
TI_liquid = U2_liquid_rms ./ abs(u2_liquid_mean);

%% Define Rose Pine colors
color_air = [0.85, 0.52, 0.60];  % rose (darker)
color_liquid = [0.56, 0.73, 0.78];  % foam

%% Plot Mean Velocity, RMS Fluctuations, and Turbulence Intensity side by side

figure;

% First subplot: Mean Velocity
subplot(1, 3, 1);
hold on;

% Plot air phase (hollow markers)
semilogx(mean_air_velocity_magnitude, y_air_line, ...
    'Color', color_air, 'Marker', 's', 'LineStyle', 'none', ...
    'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerEdgeColor', color_air, 'MarkerFaceColor', 'none', ...
    'DisplayName', 'Air Phase');

% Plot liquid phase (filled markers)
semilogx(mean_liquid_velocity_magnitude, y_liquid_line, ...
    'Color', color_liquid, 'Marker', 's', 'LineStyle', 'none', ...
    'LineWidth', 1, 'MarkerSize', 6, ...
    'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', color_liquid, ...
    'DisplayName', 'Liquid Phase');

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

% Plot air phase RMS (hollow markers)
semilogx(U1_air_rms, y_air_line, ...
    'Color', color_air, 'Marker', 's', 'LineStyle', 'none', ...
    'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerEdgeColor', color_air, 'MarkerFaceColor', 'none', ...
    'DisplayName', 'Air Phase');

% Plot liquid phase RMS (filled markers)
semilogx(U2_liquid_rms, y_liquid_line, ...
    'Color', color_liquid, 'Marker', 's', 'LineStyle', 'none', ...
    'LineWidth', 1, 'MarkerSize', 6, ...
    'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', color_liquid, ...
    'DisplayName', 'Liquid Phase');

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

% Plot air phase TI (hollow markers)
semilogx(TI_air, y_air_line, ...
    'Color', color_air, 'Marker', 's', 'LineStyle', 'none', ...
    'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerEdgeColor', color_air, 'MarkerFaceColor', 'none', ...
    'DisplayName', 'Air Phase');

% Plot liquid phase TI (filled markers)
semilogx(TI_liquid, y_liquid_line, ...
    'Color', color_liquid, 'Marker', 's', 'LineStyle', 'none', ...
    'LineWidth', 1, 'MarkerSize', 6, ...
    'MarkerEdgeColor', color_liquid, 'MarkerFaceColor', color_liquid, ...
    'DisplayName', 'Liquid Phase');

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

%% Plot velocity vector map for a single frame with center column highlighted

frame_to_plot = 1; % Choose a frame to visualize

[X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame_to_plot);

% Get the x-position of the column we're averaging
x_column = X1(1, column_index);

% Shift all x-coordinates so that the averaging column is at x = 0
X1_shifted = X1 - x_column;
X2_shifted = X2 - x_column;
X3_shifted = X3 - x_column;

figure;
hold on;

% Plot liquid phase vectors (foam color)
quiver(X2_shifted, Y2, U2, V2, 'Color', color_liquid, 'LineWidth', 1, 'AutoScale', 'on', 'AutoScaleFactor', 1.5);

% Plot air phase vectors (rose color)
quiver(X1_shifted, Y1, U1, V1, 'Color', color_air, 'LineWidth', 1, 'AutoScale', 'on', 'AutoScaleFactor', 1.5);

% Plot film height (gold/amber Rose Pine color)
color_film = [0.90, 0.77, 0.53];  % gold (same as conditional averaging)
plot(X3_shifted, Y3, '-', 'Color', color_film, 'LineWidth', 3, 'DisplayName', 'Film Height');

% Highlight the center column (now at x = 0)
xline(0, 'k--', 'LineWidth', 2, 'DisplayName', 'Averaging Column (x = 0)');

% Add labels and formatting
xlabel('x position (mm)', 'FontSize', 12);
ylabel('y position (mm)', 'FontSize', 12);
grid on;
axis equal;
ylim([0 28]);
yticks(0:2:28);
xlim([-15 15.1]);

% Add legend
legend('Liquid Phase', 'Air Phase', 'Film Height', 'Averaging Column (x = 0)', 'Location', 'best', 'FontSize', 9);
hold off;