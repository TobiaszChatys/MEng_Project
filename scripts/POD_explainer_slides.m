%% POD_EXPLAINER_SLIDES.m
%  Generates 3 presentation-quality figures explaining POD
%  Run: >> POD_explainer_slides
%  Outputs: slide1_what_is_pod.png, slide2_energy_ranking.png, slide3_filtering.png
%
%  Export tip: figures are sized 13.3" x 7.5" (standard 16:9 slide)
%              at 300 DPI — drop directly into PowerPoint.

clear; clc; close all;

%% ======================================================================
%  Colour palette (Catppuccin-inspired, high contrast for projection)
%  ======================================================================
col_bg      = [0.12 0.12 0.18];   % dark background
col_text    = [0.80 0.84 0.96];   % light text
col_dim     = [0.35 0.36 0.44];   % muted grey
col_blue    = [0.54 0.71 0.98];   % accent blue
col_green   = [0.65 0.89 0.63];   % accent green
col_red     = [0.95 0.55 0.66];   % accent red
col_yellow  = [0.98 0.89 0.68];   % accent yellow
col_mauve   = [0.80 0.65 0.97];   % accent purple

slide_w = 13.3;  % inches (16:9)
slide_h = 7.5;

%% ======================================================================
%  SLIDE 1 — "What is POD?" (equation + PCA scatter analogy)
%  ======================================================================
fig1 = figure('Color', col_bg, 'Units', 'inches', 'Position', [1 1 slide_w slide_h]);
set(fig1, 'InvertHardcopy', 'off');  % preserve dark background on export

% ----- Title -----
annotation('textbox', [0.02 0.88 0.96 0.10], ...
  'String', 'What is POD?', ...
  'FontSize', 36, 'FontWeight', 'bold', 'Color', col_text, ...
  'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');

% ----- Subtitle -----
annotation('textbox', [0.02 0.82 0.96 0.07], ...
  'String', 'Proper Orthogonal Decomposition  —  same idea as PCA, applied to flow fields', ...
  'FontSize', 18, 'Color', col_dim, ...
  'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');

% ----- Decomposition equation (left half) -----
annotation('textbox', [0.03 0.52 0.50 0.12], ...
  'String', 'Decomposition:', ...
  'FontSize', 20, 'FontWeight', 'bold', 'Color', col_blue, ...
  'EdgeColor', 'none', 'FontName', 'Helvetica');

annotation('textbox', [0.03 0.40 0.50 0.14], ...
  'String', '$$\mathbf{u}(\mathbf{x},t) = \bar{\mathbf{u}}(\mathbf{x}) + \sum_{k=1}^{N} a_k(t)\,\boldsymbol{\phi}_k(\mathbf{x})$$', ...
  'Interpreter', 'latex', ...
  'FontSize', 22, 'Color', col_text, ...
  'EdgeColor', 'none', 'FontName', 'Helvetica');

% Legend for equation terms
terms   = {'u(x,t) = original snapshot', ...
  'u-bar   = time-averaged mean', ...
  'a_k(t)  = temporal coefficients', ...
  '\phi_k(x) = spatial modes (ranked by energy)'};
colors  = {col_red, col_dim, col_yellow, col_green};
for i = 1:4
  annotation('textbox', [0.06 0.36 - i*0.065 0.44 0.06], ...
    'String', terms{i}, ...
    'FontSize', 15, 'Color', colors{i}, ...
    'EdgeColor', 'none', 'FontName', 'Helvetica');
end

% ----- PCA scatter plot analogy (right half) -----
ax1 = axes('Position', [0.55 0.08 0.42 0.68], 'Color', col_bg);
hold on; axis equal; box off;
set(ax1, 'XColor', col_dim, 'YColor', col_dim, 'FontSize', 12, 'Color', col_bg);

rng(7);
mu = [0 0];
sigma = [2 1.5; 1.5 2];
pts = mvnrnd(mu, sigma, 200);

scatter(pts(:,1), pts(:,2), 25, col_blue, 'filled', 'MarkerFaceAlpha', 0.4);

% Principal axes
[V, D] = eig(sigma);
[~, idx] = sort(diag(D), 'descend');
V = V(:, idx);

% PC1 (most energy)
quiver(0, 0, V(1,1)*3.2, V(2,1)*3.2, 0, ...
  'Color', col_green, 'LineWidth', 3.5, 'MaxHeadSize', 0.4);
text(V(1,1)*3.4, V(2,1)*3.4, {'Mode 1', '(most energy)'}, ...
  'Color', col_green, 'FontSize', 14, 'FontWeight', 'bold', ...
  'FontName', 'Helvetica');

% PC2
quiver(0, 0, V(1,2)*2.0, V(2,2)*2.0, 0, ...
  'Color', col_yellow, 'LineWidth', 3.5, 'MaxHeadSize', 0.4);
text(V(1,2)*2.2, V(2,2)*2.2, 'Mode 2', ...
  'Color', col_yellow, 'FontSize', 14, 'FontWeight', 'bold', ...
  'FontName', 'Helvetica');

xlabel('V_x', 'Color', col_dim, 'FontSize', 14);
ylabel('V_y', 'Color', col_dim, 'FontSize', 14);
title('"Rotate axes to align with maximum variance"', ...
  'Color', col_mauve, 'FontSize', 16, 'FontName', 'Helvetica');
xlim([-5 5]); ylim([-5 5]);
grid on; set(ax1, 'GridColor', col_dim, 'GridAlpha', 0.3);
hold off;

% Export
exportgraphics(fig1, 'slide1_what_is_pod.png', 'Resolution', 300);
fprintf('✓ Saved slide1_what_is_pod.png\n');


%% ======================================================================
%  SLIDE 2 — Energy ranking ("How many modes do we need?")
%  ======================================================================
fig2 = figure('Color', col_bg, 'Units', 'inches', 'Position', [1 1 slide_w slide_h]);
set(fig2, 'InvertHardcopy', 'off');

annotation('textbox', [0.02 0.88 0.96 0.10], ...
  'String', 'Modes Ranked by Energy', ...
  'FontSize', 36, 'FontWeight', 'bold', 'Color', col_text, ...
  'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');

annotation('textbox', [0.02 0.82 0.96 0.07], ...
  'String', 'First few modes capture the real physics  —  the tail is just noise', ...
  'FontSize', 18, 'Color', col_dim, ...
  'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');

% ----- Bar chart (left) -----
ax2a = axes('Position', [0.06 0.10 0.42 0.68], 'Color', col_bg);
hold on;

n_modes = 20;
% Simulated eigenvalue spectrum (realistic exponential decay)
energies = 38 * exp(-0.28 * (0:n_modes-1)) + 0.5;
energies = energies / sum(energies) * 100;  % percentage

bar_colors = zeros(n_modes, 3);
for i = 1:n_modes
  if i <= 5
    bar_colors(i,:) = col_green;   % physics
  elseif i <= 12
    bar_colors(i,:) = col_yellow;  % transition
  else
    bar_colors(i,:) = col_red;     % noise
  end
end

for i = 1:n_modes
  bar(i, energies(i), 0.7, 'FaceColor', bar_colors(i,:), ...
    'EdgeColor', col_dim, 'LineWidth', 0.5, 'FaceAlpha', 0.85);
end

% Braces / annotations
% "Real physics" label
annotation('textbox', [0.08 0.72 0.15 0.06], ...
  'String', 'real physics ↓', ...
  'FontSize', 14, 'FontWeight', 'bold', 'Color', col_green, ...
  'EdgeColor', 'none', 'FontName', 'Helvetica');

% "Noise" label
annotation('textbox', [0.34 0.72 0.10 0.06], ...
  'String', 'noise ↓', ...
  'FontSize', 14, 'FontWeight', 'bold', 'Color', col_red, ...
  'EdgeColor', 'none', 'FontName', 'Helvetica');

set(ax2a, 'XColor', col_dim, 'YColor', col_dim, 'FontSize', 12, 'Color', col_bg);
xlabel('Mode number', 'Color', col_dim, 'FontSize', 14);
ylabel('% of total energy', 'Color', col_dim, 'FontSize', 14);
title('Individual mode energy', 'Color', col_text, 'FontSize', 16, 'FontName', 'Helvetica');
xlim([0 n_modes+1]);
grid on; set(ax2a, 'GridColor', col_dim, 'GridAlpha', 0.3);
hold off;

% ----- Cumulative energy curve (right) -----
ax2b = axes('Position', [0.56 0.10 0.40 0.68], 'Color', col_bg);
hold on;

cum_energy = cumsum(energies);
plot(1:n_modes, cum_energy, '-o', 'Color', col_blue, 'LineWidth', 2.5, ...
  'MarkerFaceColor', col_blue, 'MarkerSize', 5);

% 90% threshold line
yline(90, '--', 'Color', col_mauve, 'LineWidth', 2, 'Alpha', 0.8);
% Find where it crosses 90%
idx90 = find(cum_energy >= 90, 1);
plot(idx90, cum_energy(idx90), 'o', 'MarkerSize', 14, ...
  'MarkerEdgeColor', col_mauve, 'LineWidth', 2.5);
text(idx90 + 0.8, cum_energy(idx90) - 3, ...
  sprintf('90%% energy\nat mode %d', idx90), ...
  'Color', col_mauve, 'FontSize', 14, 'FontWeight', 'bold', ...
  'FontName', 'Helvetica');

% Shade the "keep" region
fill([0 idx90 idx90 0], [0 0 100 100], col_green, ...
  'FaceAlpha', 0.08, 'EdgeColor', 'none');

% Shade the "discard" region
fill([idx90 n_modes n_modes idx90], [0 0 100 100], col_red, ...
  'FaceAlpha', 0.08, 'EdgeColor', 'none');

text(idx90/2, 15, 'KEEP', 'Color', col_green, 'FontSize', 16, ...
  'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');
text((idx90 + n_modes)/2, 15, 'DISCARD', 'Color', col_red, 'FontSize', 16, ...
  'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');

set(ax2b, 'XColor', col_dim, 'YColor', col_dim, 'FontSize', 12, 'Color', col_bg);
xlabel('Number of modes', 'Color', col_dim, 'FontSize', 14);
ylabel('Cumulative energy (%)', 'Color', col_dim, 'FontSize', 14);
title('Cumulative energy spectrum', 'Color', col_text, 'FontSize', 16, 'FontName', 'Helvetica');
xlim([0 n_modes+1]); ylim([0 105]);
grid on; set(ax2b, 'GridColor', col_dim, 'GridAlpha', 0.3);
hold off;

exportgraphics(fig2, 'slide2_energy_ranking.png', 'Resolution', 300);
fprintf('✓ Saved slide2_energy_ranking.png\n');


%% ======================================================================
%  SLIDE 3 — Before & After: Raw vs POD-filtered
%  ======================================================================
fig3 = figure('Color', col_bg, 'Units', 'inches', 'Position', [1 1 slide_w slide_h]);
set(fig3, 'InvertHardcopy', 'off');

annotation('textbox', [0.02 0.88 0.96 0.10], ...
  'String', 'POD Filtering: Before & After', ...
  'FontSize', 36, 'FontWeight', 'bold', 'Color', col_text, ...
  'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');

annotation('textbox', [0.02 0.82 0.96 0.07], ...
  'String', 'Reconstruct using top modes only  →  noise removed, physics preserved', ...
  'FontSize', 18, 'Color', col_dim, ...
  'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');

% Generate synthetic "wave-like" velocity field
rng(42);
nx = 16; ny = 10;
[X, Y] = meshgrid(linspace(0, 4, nx), linspace(0, 2, ny));

% Base flow: shear profile + sinusoidal wave
Vx_clean = 0.3 * Y + 0.15 * sin(2*pi*X/2.5) .* exp(-0.5*(Y-1).^2);
Vy_clean = 0.08 * cos(2*pi*X/2.5) .* exp(-0.5*(Y-1).^2);

% Wave interface (sinusoidal surface)
x_interface = linspace(0, 4, 100);
y_interface = 1.3 + 0.2 * sin(2*pi*x_interface/2.5);

% ----- LEFT: Raw (noisy) -----
ax3a = axes('Position', [0.04 0.08 0.40 0.70], 'Color', col_bg);
hold on;

Vx_noisy = Vx_clean + 0.15 * randn(ny, nx);
Vy_noisy = Vy_clean + 0.15 * randn(ny, nx);

% Add extra-bad vectors near interface (y ≈ 1.0–1.5)
interface_mask = (Y > 0.8) & (Y < 1.6);
Vx_noisy(interface_mask) = Vx_noisy(interface_mask) + 0.3 * randn(sum(interface_mask(:)), 1);
Vy_noisy(interface_mask) = Vy_noisy(interface_mask) + 0.3 * randn(sum(interface_mask(:)), 1);

% Velocity magnitude for colour
mag_noisy = sqrt(Vx_noisy.^2 + Vy_noisy.^2);
quiver(X, Y, Vx_noisy, Vy_noisy, 0.8, 'Color', col_red, 'LineWidth', 1.2);

% Interface line
plot(x_interface, y_interface, '--', 'Color', col_yellow, 'LineWidth', 2);
text(3.5, 1.7, 'interface', 'Color', col_yellow, 'FontSize', 11, 'FontName', 'Helvetica');

% Highlight bad vectors
annotation('textbox', [0.10 0.45 0.20 0.06], ...
  'String', '← spurious vectors', ...
  'FontSize', 12, 'Color', col_red, 'FontWeight', 'bold', ...
  'EdgeColor', 'none', 'FontName', 'Helvetica');

set(ax3a, 'XColor', col_dim, 'YColor', col_dim, 'FontSize', 11, 'Color', col_bg);
xlabel('x (mm)', 'Color', col_dim, 'FontSize', 12);
ylabel('y (mm)', 'Color', col_dim, 'FontSize', 12);
title('Raw PIV snapshot', 'Color', col_red, 'FontSize', 18, ...
  'FontWeight', 'bold', 'FontName', 'Helvetica');
xlim([-0.3 4.3]); ylim([-0.3 2.3]);
grid on; set(ax3a, 'GridColor', col_dim, 'GridAlpha', 0.2);
hold off;

% ----- ARROW between panels -----
annotation('textarrow', [0.46 0.52], [0.45 0.45], ...
  'Color', col_yellow, 'LineWidth', 3, 'HeadWidth', 15, 'HeadLength', 12);
annotation('textbox', [0.44 0.48 0.12 0.06], ...
  'String', 'POD filter', ...
  'FontSize', 14, 'FontWeight', 'bold', 'Color', col_yellow, ...
  'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');
annotation('textbox', [0.44 0.38 0.12 0.06], ...
  'String', '(90% energy)', ...
  'FontSize', 11, 'Color', col_dim, ...
  'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
  'FontName', 'Helvetica');

% ----- RIGHT: POD-filtered (clean) -----
ax3b = axes('Position', [0.56 0.08 0.40 0.70], 'Color', col_bg);
hold on;

% Add very slight residual noise to look realistic (not perfectly smooth)
Vx_filtered = Vx_clean + 0.015 * randn(ny, nx);
Vy_filtered = Vy_clean + 0.015 * randn(ny, nx);

quiver(X, Y, Vx_filtered, Vy_filtered, 0.8, 'Color', col_green, 'LineWidth', 1.2);

% Interface line
plot(x_interface, y_interface, '--', 'Color', col_yellow, 'LineWidth', 2);

% Annotations
annotation('textbox', [0.62 0.45 0.25 0.06], ...
  'String', 'coherent structures preserved →', ...
  'FontSize', 12, 'Color', col_green, 'FontWeight', 'bold', ...
  'EdgeColor', 'none', 'FontName', 'Helvetica');

set(ax3b, 'XColor', col_dim, 'YColor', col_dim, 'FontSize', 11, 'Color', col_bg);
xlabel('x (mm)', 'Color', col_dim, 'FontSize', 12);
ylabel('y (mm)', 'Color', col_dim, 'FontSize', 12);
title('POD-filtered (90% energy)', 'Color', col_green, 'FontSize', 18, ...
  'FontWeight', 'bold', 'FontName', 'Helvetica');
xlim([-0.3 4.3]); ylim([-0.3 2.3]);
grid on; set(ax3b, 'GridColor', col_dim, 'GridAlpha', 0.2);
hold off;

exportgraphics(fig3, 'slide3_filtering.png', 'Resolution', 300);
fprintf('✓ Saved slide3_filtering.png\n');

fprintf('\n=== All 3 slides exported as 16:9 PNGs at 300 DPI ===\n');
fprintf('Drop them straight into PowerPoint as full-slide images.\n');
