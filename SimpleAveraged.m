clc; clear; close all;

[S, filename] = loadData('L8_G9.mat');

frames = size(S.smoothed_film_height_matrix_out, 2) - 1;


% picking a fixed column to perform simple time-averaging on
column_index = 25;

% limit data to valid index:
U1_air_line = squeeze(S.all_u_matrix_air(:, column_index, 1:frames));
U2_liquid_line = squeeze(S.all_u_matrix_liquid(:, column_index, 1:frames));

u1_air_mean = mean(U1_air_line, 1, 'omitnan');
u2_liquid_mean = mean(U2_liquid_line, 1, 'omitnan');



