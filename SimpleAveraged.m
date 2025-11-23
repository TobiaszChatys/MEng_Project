clc; clear; close all;

[S, filename] = loadData('L8_G9.mat');

frames = size(S.smoothed_film_height_matrix_out, 2) - 1;


% picking a fixed column to perform simple time-averaging on
column_index = 25;

% collect data across all frames:
U1_all = S.all_u_matrix_air(:, column_index, 1:frames);



