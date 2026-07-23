%% slide8_mode_analysis.m
%  Produces Slide 8 — Mode Pair Phase Plot (presentation-ready)
%
%  Single figure (16:9): temporal coefficients for a mode pair overlaid,
%  showing the ~90° phase shift that proves travelling-wave structure.
%  Phase shift computed via Hilbert transform, annotated on the plot.
%
%  Run from MEng_Project root:
%    >> run scripts/slide8_mode_analysis.m
%
%  Output:  scripts/slide8_output/<CASE>_mode_pair.png
%           scripts/slide8_output/<CASE>_mode_pair.fig

tic; clc; close all;

%% ── Settings ─────────────────────────────────────────────────────────
CASE_NAME   = 'L8_G6';
MODE_A      = 1;               % first mode in pair
MODE_B      = 2;               % second mode in pair
XLIM_SEC    = [0 0.5];        % time window to display (seconds)
FIG_SIZE    = [50 50 1920 1080];

SCRIPT_DIR  = fileparts(mfilename('fullpath'));
PROJ_ROOT   = fileparts(SCRIPT_DIR);
addpath(fullfile(PROJ_ROOT, 'src'));

% sampling rate comes from the case config, not all cases are 1000hz
case_info   = caseConfig(CASE_NAME);
FS          = case_info.fs;
DATA_DIR    = fullfile(PROJ_ROOT, 'data', 'Cases');
POD_DIR     = fullfile(SCRIPT_DIR, 'Results', 'POD_data');
OUT_DIR     = fullfile(SCRIPT_DIR, 'slide8_output');
if ~exist(POD_DIR, 'dir'), mkdir(POD_DIR); end
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 1 — Load POD results
%% ═══════════════════════════════════════════════════════════════════════
pod_file = fullfile(POD_DIR, ['POD_Results_' CASE_NAME '.mat']);
if ~exist(pod_file, 'file')
  error('POD results not found at %s. Run slide5 or Incremental_POD first.', pod_file);
end

fprintf('Loading POD results for %s...\n', CASE_NAME);
d = load(pod_file);

% Eigenvectors (temporal coefficients)
if isfield(d, 'eigenvectors_matrix'), eigvecs = double(d.eigenvectors_matrix);
elseif isfield(d, 'eigvecs'), eigvecs = double(d.eigvecs);
else, error('Eigenvectors not found in POD file'); end

if isfield(d, 'sort_index'), sort_idx = double(d.sort_index(:));
elseif isfield(d, 'sort_idx'), sort_idx = double(d.sort_idx(:));
else, error('Sort index not found in POD file'); end

if isfield(d, 'eigenvalues'), eigenvalues = double(d.eigenvalues(:));
elseif isfield(d, 'eigvals_sorted'), eigenvalues = double(d.eigvals_sorted(:));
else, eigenvalues = []; end

clear d

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 2 — Extract temporal coefficients & compute phase shift
%% ═══════════════════════════════════════════════════════════════════════
N_frames = size(eigvecs, 1);
time_vec = (1:N_frames) / FS;   % time in seconds

a_A = eigvecs(:, sort_idx(MODE_A));
a_B = eigvecs(:, sort_idx(MODE_B));

% Energy percentages
if ~isempty(eigenvalues)
  total_E = sum(eigenvalues);
  pct_A = eigenvalues(MODE_A) / total_E * 100;
  pct_B = eigenvalues(MODE_B) / total_E * 100;
else
  pct_A = NaN; pct_B = NaN;
end

% Hilbert phase shift (FFT-based, no toolbox)
H_A = hilbert_fft(a_A);
H_B = hilbert_fft(a_B);
phase_diff = rad2deg(angle(H_A ./ H_B));
median_phase = median(phase_diff);
fprintf('  Mode %d vs %d: median phase shift = %.1f°\n', MODE_A, MODE_B, median_phase);

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 3 — Presentation-ready figure
%% ═══════════════════════════════════════════════════════════════════════
TEXT_COL   = [0.13 0.13 0.13];
col_A      = [0.0 0.45 0.85];   % blue
col_B      = [0.0 0.75 0.45];   % green
FONT_LABEL = 22;
FONT_TITLE = 26;
FONT_TICK  = 16;
FONT_LEG   = 18;

fig = figure('Color', 'w', 'Position', FIG_SIZE, 'Visible', 'off');
ax  = axes(fig);
hold(ax, 'on');

plot(ax, time_vec, a_A, '-', 'Color', col_A, 'LineWidth', 2.2, ...
  'DisplayName', sprintf('Mode %d  (%.1f\\%% energy)', MODE_A, pct_A));
plot(ax, time_vec, a_B, '-', 'Color', col_B, 'LineWidth', 2.2, ...
  'DisplayName', sprintf('Mode %d  (%.1f\\%% energy)', MODE_B, pct_B));

hold(ax, 'off');

set(ax, ...
  'XLim',                XLIM_SEC, ...
  'FontSize',            FONT_TICK, ...
  'FontName',            'Helvetica', ...
  'TickDir',             'in', ...
  'TickLabelInterpreter','latex', ...
  'LineWidth',           1.4, ...
  'Box',                 'on', ...
  'XColor',              TEXT_COL, ...
  'YColor',              TEXT_COL);

xlabel(ax, 'Time (s)', 'Interpreter', 'latex', 'FontSize', FONT_LABEL);
ylabel(ax, 'Temporal Coefficient', 'Interpreter', 'latex', 'FontSize', FONT_LABEL);
title(ax, sprintf('Mode %d vs Mode %d %s', MODE_A, MODE_B, ...
  strrep(CASE_NAME, '_', '_')), ...
  'Interpreter', 'latex', 'FontSize', FONT_TITLE, 'Color', TEXT_COL);

grid(ax, 'on');
ax.GridAlpha = 0.15;
ax.GridColor = [0.5 0.5 0.5];

% Legend
leg = legend(ax, 'Location', 'northeast', 'Interpreter', 'latex', ...
  'FontSize', FONT_LEG, 'Box', 'on');
leg.EdgeColor = [0.6 0.6 0.6];

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 4 — Export
%% ═══════════════════════════════════════════════════════════════════════
fig_path = fullfile(OUT_DIR, [CASE_NAME '_mode_pair.fig']);
saveas(fig, fig_path);
fprintf('Saved:  %s\n', fig_path);

png_path = fullfile(OUT_DIR, [CASE_NAME '_mode_pair.png']);
exportgraphics(fig, png_path, 'Resolution', 300, 'BackgroundColor', 'white');
fprintf('Saved:  %s\n', png_path);

close(fig);
fprintf('Total elapsed: %.1f s\n', toc);

%% ── Local helper: Hilbert transform via FFT (no toolbox) ────────────
function h = hilbert_fft(x)
N = length(x);
Xf = fft(x);
H = zeros(N, 1);
H(1) = 1;
if mod(N,2) == 0
  H(N/2+1) = 1;
  H(2:N/2) = 2;
else
  H(2:(N+1)/2) = 2;
end
h = ifft(Xf .* H);
end
