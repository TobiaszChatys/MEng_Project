%% slide4_figures.m
%  Generates three publication-ready PNGs for Slide 4 (POD explanation).
%
%  Figures produced:
%    1. pca_analogy.png        — scatter plot with rotated principal axes
%    2. snapshot_matrix.png    — PIV grid → column vector → matrix X
%    3. eigenvalue_spectrum.png — real eigenvalue bars + cumulative energy
%
%  Run from the MEng_Project root:
%    >> run scripts/slide4_figures.m
%
%  Outputs saved to:  scripts/slide4_output/

clc; clear; close all;

%% ── Paths ────────────────────────────────────────────────────────────
project_root = fileparts(mfilename('fullpath'));
pod_file     = fullfile(project_root, 'Results', 'POD_data', 'POD_Results_L8_G12.mat');
out_dir      = fullfile(project_root, 'slide4_output');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ── Load real POD data ───────────────────────────────────────────────
fprintf('Loading POD results...\n');
data           = load(pod_file, 'eigenvalues', 'cumulative_energy', ...
  'rows_liquid',  'columns_liquid',   ...
  'number_of_modes_at_99');
eigenvalues    = double(data.eigenvalues(:));
cum_energy     = double(data.cumulative_energy(:));
n_modes_99     = double(data.number_of_modes_at_99);
grid_rows      = double(data.rows_liquid);
grid_cols      = double(data.columns_liquid);
total_modes    = numel(eigenvalues);

%% ── Colour palette (MATLAB winter colourmap) ─────────────────────────
SIGNAL  = [  0,  85, 212] / 255;   % strong blue
TEAL    = [  0, 128, 128] / 255;   % mid-winter
SEAFOAM = [  0, 179, 102] / 255;   % end-winter
NOISE   = [204,  34,  34] / 255;   % red
ENERGY  = [212, 138,   0] / 255;   % amber
MUTED   = [136, 136, 136] / 255;   % grey
TEXT    = [ 34,  34,  34] / 255;   % near-black

%% ════════════════════════════════════════════════════════════════════
%  FIGURE 1 — PCA scatter-plot analogy
%% ════════════════════════════════════════════════════════════════════
fprintf('Figure 1: PCA analogy...\n');

rng(42);
n_pts = 250;
ang   = pi / 6;
R     = [cos(ang), -sin(ang); sin(ang), cos(ang)];
raw   = [2.0 * randn(n_pts, 1), 0.5 * randn(n_pts, 1)];
pts   = (R * raw')';

fig1 = figure('Color', 'white', 'Position', [100 100 700 520]);
ax   = axes(fig1);

scatter(pts(:,1), pts(:,2), 18, TEAL, 'filled', ...
  'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
hold on;

% Dashed original-axes grid lines
xline(0, '--', 'Color', MUTED, 'LineWidth', 0.8);
yline(0, '--', 'Color', MUTED, 'LineWidth', 0.8);

% Mode 1 arrow (maximum energy direction)
m1 = R * [3.5; 0];
m2 = R * [0;   1.6];

quiver(0, 0, m1(1), m1(2), 0, 'Color', SIGNAL, ...
  'LineWidth', 2.5, 'MaxHeadSize', 0.4);
quiver(0, 0, m2(1), m2(2), 0, 'Color', NOISE,  ...
  'LineWidth', 2.5, 'MaxHeadSize', 0.6);

text(m1(1) + 0.15, m1(2) + 0.28, ...
  '$\phi_1$', ...
  'Interpreter', 'latex', 'FontSize', 16, ...
  'Color', SIGNAL, 'FontWeight', 'bold');
text(m2(1) - 0.2, m2(2) + 0.25, ...
  '$\phi_2$', ...
  'Interpreter', 'latex', 'FontSize', 16, ...
  'Color', NOISE, 'FontWeight', 'bold');

xlim([-5 5]); ylim([-3.5 4]);
axis equal;
xlabel('$V_x$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$V_y$', 'Interpreter', 'latex', 'FontSize', 18);
title('POD: rotating axes to align with maximum energy', ...
  'Interpreter', 'latex', 'FontSize', 16);

% Equation box bottom-right
eq_str = ['$u(\mathbf{x},t) = \bar{u}(\mathbf{x})' ...
  ' + \sum_{k} a_k(t)\,\phi_k(\mathbf{x})$'];
text(4.7, -3.1, eq_str, ...
  'Interpreter', 'latex', 'FontSize', 14, ...
  'HorizontalAlignment', 'right', ...
  'BackgroundColor', 'white', 'EdgeColor', MUTED, ...
  'Margin', 4, 'Color', TEXT);

style_ax(ax);

exportgraphics(fig1, fullfile(out_dir, 'pca_analogy.png'), ...
  'Resolution', 300, 'BackgroundColor', 'white');
fprintf('  saved pca_analogy.png\n');

%% ════════════════════════════════════════════════════════════════════
%  FIGURE 2 — Snapshot matrix diagram
%% ════════════════════════════════════════════════════════════════════
fprintf('Figure 2: Snapshot matrix...\n');

fig2 = figure('Color', 'white', 'Position', [100 100 1100 450]);
ax2  = axes(fig2, 'Position', [0.04 0.06 0.92 0.88]);
hold on; axis off;
xlim([-0.5 11]); ylim([-1.5 3.8]);

% Common vertical centre for all three elements
Y_MID = 1.1;

% ── PIV arrow grid (left) ─────────────────────────────────────────
rng(7);
gx = 4; gy = 4;
sp = 0.55;
grid_w = (gx - 1) * sp;
grid_h = (gy - 1) * sp;
piv_ox = 0.3;
piv_oy = Y_MID - grid_h / 2;

for i = 1:gy
  for j = 1:gx
    x0 = piv_ox + (j-1) * sp;
    y0 = piv_oy + (gy - i) * sp;
    dx = 0.22 + 0.08 * randn();
    dy = 0.07 * randn();
    quiver(x0, y0, dx, dy, 0, 'Color', SIGNAL, ...
      'LineWidth', 1.2, 'MaxHeadSize', 0.6, 'AutoScale', 'off');
  end
end

% Rectangle border around grid
pad = 0.25;
piv_rx = piv_ox - pad;
piv_ry = piv_oy - pad;
piv_rw = grid_w + 2 * pad;
piv_rh = grid_h + 2 * pad;
rectangle('Position', [piv_rx piv_ry piv_rw piv_rh], 'EdgeColor', MUTED, ...
  'LineWidth', 1.2, 'Curvature', [0.05 0.05]);

piv_cx = piv_ox + grid_w / 2;
text(piv_cx, piv_ry - 0.25, 'PIV snapshot', ...
  'FontSize', 12, 'Color', MUTED, 'HorizontalAlignment', 'center', ...
  'Interpreter', 'none');
text(piv_cx, piv_ry - 0.55, ...
  sprintf('(%d $\\times$ %d grid)', grid_rows, grid_cols), ...
  'FontSize', 11, 'Color', MUTED, 'HorizontalAlignment', 'center', ...
  'Interpreter', 'latex');

% ── Arrow: PIV → column vector ────────────────────────────────────
arrow_y = Y_MID;
arr1_x1 = piv_rx + piv_rw + 0.15;
arr1_x2 = arr1_x1 + 0.7;
n1 = data2norm(ax2, [arr1_x1 arr1_x2], [arrow_y arrow_y]);
annotation(fig2, 'arrow', n1(1,:), n1(2,:), ...
  'Color', TEXT, 'HeadWidth', 8, 'HeadLength', 8, 'LineWidth', 1.5);

% ── Column vector ─────────────────────────────────────────────────
entries = {'$V_{x_1}$', '$V_{y_1}$', '$\vdots$', '$V_{x_n}$', '$V_{y_n}$'};
n_e     = numel(entries);
row_h   = 0.42;
col_h   = n_e * row_h;
col_x   = arr1_x2 + 0.3;
col_bot = Y_MID - col_h / 2;

rectangle('Position', [col_x - 0.1, col_bot - 0.06, 0.75, col_h + 0.12], ...
  'EdgeColor', TEXT, 'FaceColor', [0.94 0.97 1.0], ...
  'LineWidth', 1, 'Curvature', [0.05 0.05]);

for k = 1:n_e
  text(col_x + 0.27, col_bot + (n_e - k + 0.5) * row_h, entries{k}, ...
    'Interpreter', 'latex', 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

text(col_x + 0.27, col_bot - 0.35, '$\mathbf{x}_k$', ...
  'Interpreter', 'latex', 'FontSize', 14, ...
  'HorizontalAlignment', 'center', 'Color', SIGNAL, 'FontWeight', 'bold');

% ── Arrow: column vector → matrix ─────────────────────────────────
arr2_x1 = col_x + 0.75 + 0.15;
arr2_x2 = arr2_x1 + 0.7;
n2 = data2norm(ax2, [arr2_x1 arr2_x2], [arrow_y arrow_y]);
annotation(fig2, 'arrow', n2(1,:), n2(2,:), ...
  'Color', TEXT, 'HeadWidth', 8, 'HeadLength', 8, 'LineWidth', 1.5);

% ── Matrix X ──────────────────────────────────────────────────────
mat_x   = arr2_x2 + 0.6;
mat_w   = 3.8;
mat_h   = 2.6;
mat_bot = Y_MID - mat_h / 2;

% "X =" label
text(mat_x - 0.35, Y_MID, '$X =$', ...
  'Interpreter', 'latex', 'FontSize', 18, ...
  'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
  'Color', ENERGY, 'FontWeight', 'bold');

rectangle('Position', [mat_x, mat_bot, mat_w, mat_h], ...
  'EdgeColor', TEXT, 'FaceColor', [0.94 0.97 1.0], ...
  'LineWidth', 1.2, 'Curvature', [0.03 0.05]);

% Column dividers + labels (evenly spaced inside matrix)
col_labels = {'$\mathbf{x}_1$', '$\mathbf{x}_2$', '$\cdots$', '$\mathbf{x}_{2500}$'};
n_col      = numel(col_labels);
col_gap    = mat_w / n_col;

for k = 1:n_col
  % Label centred within its column
  lx = mat_x + (k - 0.5) * col_gap;
  text(lx, Y_MID, col_labels{k}, ...
    'Interpreter', 'latex', 'FontSize', 13, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
  % Divider line on the right edge of the column (skip last)
  if k < n_col
    dx = mat_x + k * col_gap;
    line([dx dx], [mat_bot + 0.3, mat_bot + mat_h - 0.3], ...
      'Color', [SIGNAL, 0.35], 'LineWidth', 0.8);
  end
end

% Double-arrow below matrix
text(mat_x + mat_w / 2, mat_bot - 0.45, '2,500 snapshots', ...
  'FontSize', 11, 'Color', MUTED, 'HorizontalAlignment', 'center');
nda = data2norm(ax2, [mat_x + 0.1, mat_x + mat_w - 0.1], ...
  [mat_bot - 0.15, mat_bot - 0.15]);
annotation(fig2, 'doublearrow', nda(1,:), nda(2,:), ...
  'Color', MUTED, 'LineWidth', 0.8);

% Dimension note on right
dim_label = sprintf('$%d\\,\\times\\,2500$', 2 * grid_rows * grid_cols);
text(mat_x + mat_w + 0.2, Y_MID, dim_label, ...
  'Interpreter', 'latex', 'FontSize', 11, 'Color', MUTED, ...
  'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
  'Rotation', 90);

exportgraphics(fig2, fullfile(out_dir, 'snapshot_matrix.png'), ...
  'Resolution', 300, 'BackgroundColor', 'white');
fprintf('  saved snapshot_matrix.png\n');

%% ════════════════════════════════════════════════════════════════════
%  FIGURE 3 — Eigenvalue spectrum (real data)
%% ════════════════════════════════════════════════════════════════════
fprintf('Figure 3: Eigenvalue spectrum...\n');

N     = 50;
eig_n = eigenvalues(1:N) ./ eigenvalues(1);   % normalised

% Winter colourmap for bars: blue → teal → seafoam
cmap = zeros(N, 3);
for k = 1:N
  t = (k - 1) / (N - 1);
  if t < 0.5
    cmap(k,:) = (1 - 2*t) * SIGNAL  + 2*t * TEAL;
  else
    cmap(k,:) = (1 - 2*(t-0.5)) * TEAL + 2*(t-0.5) * SEAFOAM;
  end
end

fig3 = figure('Color', 'white', 'Position', [100 100 900 500]);
ax3  = axes(fig3);

for k = 1:N
  bar(k, eig_n(k), 'FaceColor', cmap(k,:), 'EdgeColor', 'white', ...
    'LineWidth', 0.3, 'BarWidth', 0.8);
  hold on;
end
ylim([0 1.12]);
ylabel('Normalised eigenvalue $(\lambda_k / \lambda_1)$', ...
  'Interpreter', 'latex', 'FontSize', 16, 'Color', TEXT);
xlim([0 N+1]);
xticks([1, 5:5:N]);
xlabel('Mode number', 'Interpreter', 'latex', 'FontSize', 16);
style_ax(ax3);

% Energy ranking equation — top-right corner
text(N - 1, 1.06, '$\lambda_1 > \lambda_2 > \cdots > \lambda_n$', ...
  'Interpreter', 'latex', 'FontSize', 16, 'Color', TEXT, ...
  'HorizontalAlignment', 'right', ...
  'BackgroundColor', 'white', 'EdgeColor', MUTED, 'Margin', 4);


exportgraphics(fig3, fullfile(out_dir, 'eigenvalue_spectrum.png'), ...
  'Resolution', 300, 'BackgroundColor', 'white');
fprintf('  saved eigenvalue_spectrum.png\n');

fprintf('\nAll figures saved to:  %s\n', out_dir);

%% ── Local functions ──────────────────────────────────────────────────
function style_ax(ax)
set(ax, ...
  'FontSize',        14, ...
  'FontName',        'Helvetica', ...
  'TickLabelInterpreter', 'latex', ...
  'XMinorTick',      'off', ...
  'YMinorTick',      'off', ...
  'LineWidth',       1.2, ...
  'XColor',          [34 34 34]/255, ...
  'YColor',          [34 34 34]/255);
end

function nxy = data2norm(ax, xd, yd)
% Convert data coordinates to normalised figure coordinates
% Returns 2×N matrix: row 1 = x_norm, row 2 = y_norm
pos = ax.Position;
xl  = ax.XLim;
yl  = ax.YLim;
nx  = pos(1) + (xd - xl(1)) / (xl(2) - xl(1)) * pos(3);
ny  = pos(2) + (yd - yl(1)) / (yl(2) - yl(1)) * pos(4);
nxy = [nx; ny];
end
