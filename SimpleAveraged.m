clc; clear; close all;

[S, filename] = loadData('L8_G9.mat');

frames = size(S.smoothed_film_height_matrix_out, 2) - 1;


% picking a fixed column to perform simple time-averaging on
column_index = 25;

% limit data to valid index:
U1_air_line = squeeze(S.all_u_matrix_air(:, column_index, 1:frames));
U2_liquid_line = squeeze(S.all_u_matrix_liquid(:, column_index, 1:frames));
V1_air_line = squeeze(S.all_v_matrix_air(:, column_index, 1:frames));
V2_liquid_line = squeeze(S.all_v_matrix_liquid(:, column_index, 1:frames));



u1_air_mean = mean(U1_air_line, 1, 'omitnan');
v1_air_mean = mean(V1_air_line, 1, 'omitnan');
u2_liquid_mean = mean(U2_liquid_line, 1, 'omitnan');
v2_liquid_mean = mean(V2_liquid_line, 1, 'omitnan');


mean_air_velocity_magnitude = hypot(u1_air_mean, v1_air_mean);
mean_liquid_velocity_magnitude = hypot(u2_liquid_mean, v2_liquid_mean);


