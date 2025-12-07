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

% jplotting the simple time-averaged velocity profiles
figure;
plot(mean_air_velocity_magnitude, y_air_line, 'r-', 'LineWidth', 2);
hold on;
plot(mean_liquid_velocity_magnitude, y_liquid_line, 'b-', 'LineWidth', 2);
ylabel('Y Position (mm)');
xlabel('Mean Velocity Magnitude');
title('Simple Time-Averaged Velocity Profiles at X = 0 mm');
legend('Air Phase', 'Liquid Phase');
grid on;
