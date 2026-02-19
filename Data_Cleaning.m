clc; clear; close all;

%% Load Data
frame = 1; % Example frame index
[S, filename] = loadData('L8_G9.mat');

X2 = S.all_transposed_x_position_matrix_liquid(:, :, frame) * 1e3;
Liquid_height = S.all_transposed_y_position_matrix_liquid(:, :, frame) * 1e3;
U2 = S.all_u_matrix_liquid(:, :, frame);
V2 = S.all_v_matrix_liquid(:, :, frame);

% Set zero velocities to NaN
zero_mask = (U2 == 0) & (V2 == 0);
U2(zero_mask) = NaN;
V2(zero_mask) = NaN;

Z2 = hypot(U2, V2);

X3 = S.film_x_position(1, :) * 1e3;
Y3 = S.smoothed_film_height_matrix_out(:, frame) * 1e3;

remove_negative_indices = Y3 >= 0;
X3 = X3(remove_negative_indices);
Y3 = Y3(remove_negative_indices);

% For any liquid vertical positon (Y2) that is less than the film height at the same x positon,
% set the liquid velocity to NaN

x_liquid = X2(1, :); % 1 x 50 (x positions)
film_height_at_liquid_x = interp1(X3, Y3, x_liquid, 'linear', 'extrap'); % 1 x 50 

liquid_above_interface = Liquid_height > film_height_at_liquid_x; % 99 x 50 logical matrix

New_U2 = U2;
New_U2(liquid_above_interface) = NaN;
V2(liquid_above_interface) = NaN;


count = sum(liquid_above_interface, 'all');