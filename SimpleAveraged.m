clc; clear; close all;

[S, filename] = loadData('L8_G9.mat');
frames = size(S.smoothed_film_height_matrix_out, 2) - 1;


% initalise storage for accumulated data

% air phase

U1_all = [];
V1_all = [];
X1_all = [];
Y1_all = [];

% liquid phase

U2_all = [];
V2_all = [];
X2_all = [];
Y2_all = [];

% film interface

X3_all = [];
Y3_all = [];

for frame = 1:frames

    [U1, V1, X1, Y1, U2, V2, X2, Y2, ~, ~, X3, Y3] = getData(S, frame);

    U1_all = [U1_all; U1(:)];
    V1_all = [V1_all; V1(:)];
    X1_all = [X1_all; X1(:)];
    Y1_all = [Y1_all; Y1(:)];

    U2_all = [U2_all; U2(:)];
    V2_all = [V2_all; V2(:)];
    X2_all = [X2_all; X2(:)];
    Y2_all = [Y2_all; Y2(:)];

    X3_all = [X3_all; X3(:)];
    Y3_all = [Y3_all; Y3(:)];
end

