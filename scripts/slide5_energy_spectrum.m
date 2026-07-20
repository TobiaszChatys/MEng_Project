%% slide5_energy_spectrum.m
%  Computes incremental POD for cases L8G2, L8G3, L8G6, L8G9 (liquid phase
%  only) and produces the Slide 5 figure:
%
%    Cumulative Energy (%) vs Mode Number — all cases overlaid
%    with 80% and 90% threshold lines, log-log scaling option.
%
%  Run from the MEng_Project root:
%    >> run scripts/slide5_energy_spectrum.m
%
%  POD results for each case are cached in scripts/Results/POD_data/ so
%  subsequent runs skip the expensive computation.
%
%  Output:  scripts/slide4_output/energy_spectrum.png

tic
clc; close all;

%% ── Settings ─────────────────────────────────────────────────────────
CASES = {
  'L8_G2',  '2 m s$^{-1}$  (2D)';
  'L8_G3',  '3 m s$^{-1}$  (Transition)';
  'L8_G6',  '6 m s$^{-1}$  (3D)';
  'L8_G9',  '9 m s$^{-1}$  (Disturbance)';
  };

SCRIPT_DIR  = fileparts(mfilename('fullpath'));
PROJ_ROOT   = fileparts(SCRIPT_DIR);
addpath(genpath(fullfile(PROJ_ROOT, 'src')));

DATA_DIR    = fullfile(PROJ_ROOT,  'data', 'Cases');
RESULTS_DIR = fullfile(SCRIPT_DIR, 'Results', 'POD_data');
OUT_DIR     = fullfile(SCRIPT_DIR, 'slide4_output');
if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end
if ~exist(OUT_DIR,     'dir'), mkdir(OUT_DIR);     end

BLOCK_SIZE  = 350;      % frames per block (matches Incremental_POD.m)
THRESHOLDS  = [0.80, 0.90];

%% ── Colour palette ──────────────────────────────────────────────────
% Winter colourmap: sample 4 evenly-spaced colours
WIN_T   = linspace(0, 1, size(CASES,1));
COLOURS = zeros(size(CASES,1), 3);
for k = 1:size(CASES,1)
  t = WIN_T(k);
  COLOURS(k,:) = [0, t, 1 - 0.5*t];   % R=0, G=t, B=1→0.5
end
% Override G3 with a bold highlight colour so its curve stands out
G3_IDX        = 2;
COLOURS(G3_IDX,:) = [0.85, 0.20, 0.00];   % burnt orange — clearly distinct

MUTED  = [0.55 0.55 0.55];
TEXT   = [0.13 0.13 0.13];
LINE_W = [1.8, 3.2, 1.8, 1.8];   % G3 thicker

%% ── Main loop: compute or load POD for each case ────────────────────
cum_energy_all = cell(size(CASES,1), 1);

for c = 1:size(CASES,1)
  
  case_name  = CASES{c,1};
  mat_file   = fullfile(DATA_DIR,    [case_name '.mat']);
  save_file  = fullfile(RESULTS_DIR, ['POD_Results_' case_name '.mat']);
  
  % ── Load from cache if available ─────────────────────────────────
  if exist(save_file, 'file')
    fprintf('[%s]  Loading cached POD results...\n', case_name);
    d = load(save_file);
    if isfield(d, 'cumulative_energy'), cum_energy_all{c} = double(d.cumulative_energy(:));
    elseif isfield(d, 'cum_energy'), cum_energy_all{c} = double(d.cum_energy(:));
    else, error('No cumulative energy found in %s', save_file); end
    continue
  end
  
  % ── Otherwise run incremental POD ────────────────────────────────
  fprintf('[%s]  Running incremental POD (liquid phase)...\n', case_name);
  
  % Load raw data
  tmp = load(mat_file);
  S   = tmp.S2P_PIV_full_mat_vars;
  
  frames = size(S.all_u_matrix_liquid, 3);
  [rows_liq, cols_liq, ~] = size(S.all_u_matrix_liquid);
  sp_liq = rows_liq * cols_liq;
  
  % Vectorise
  fprintf('  Vectorising %d frames...\n', frames);
  U_vec = zeros(sp_liq, frames);
  V_vec = zeros(sp_liq, frames);
  parfor fr = 1:frames
    U_vec(:,fr) = reshape(S.all_u_matrix_liquid(:,:,fr), [], 1);
    V_vec(:,fr) = reshape(S.all_v_matrix_liquid(:,:,fr), [], 1);
  end
  Snap = [U_vec; V_vec];
  Snap(isnan(Snap)) = 0;
  
  % Demean
  mean_vec  = mean(Snap, 2);
  Snap_fluc = Snap - mean_vec;
  clear U_vec V_vec Snap
  
  % Segment into blocks
  n_blocks = ceil(frames / BLOCK_SIZE);
  blocks   = cell(n_blocks, 1);
  parfor b = 1:n_blocks
    s = (b-1)*BLOCK_SIZE + 1;
    e = min(b*BLOCK_SIZE, frames);
    blocks{b} = Snap_fluc(:, s:e);
  end
  
  % Build temporal covariance incrementally
  fprintf('  Building temporal covariance (%d blocks)...\n', n_blocks);
  C = zeros(frames, frames);
  for b = 1:n_blocks
    sb = (b-1)*BLOCK_SIZE + 1;
    eb = min(b*BLOCK_SIZE, frames);
    for ub = b:n_blocks
      sub = (ub-1)*BLOCK_SIZE + 1;
      eub = min(ub*BLOCK_SIZE, frames);
      blk = blocks{b}' * blocks{ub};
      C(sb:eb, sub:eub) = blk;
      if b ~= ub
        C(sub:eub, sb:eb) = blk';
      end
    end
  end
  C = C / frames;
  clear blocks Snap_fluc
  
  % Eigenvalue decomposition
  fprintf('  Eigenvalue decomposition...\n');
  [eigvecs, eigval_mat] = eig(C);
  eigvals = diag(eigval_mat);
  [eigvals_sorted, sort_idx] = sort(eigvals, 'descend');
  eigvals_sorted = max(eigvals_sorted, 0);   % numerical floor
  
  total_energy   = sum(eigvals_sorted);
  cum_energy     = cumsum(eigvals_sorted) / total_energy;
  n99            = find(cum_energy >= 0.99, 1);
  
  cum_energy_all{c} = cum_energy;
  
  % Save — eigenvalues only (no spatial modes to keep file small)
  fprintf('  Saving results to %s\n', save_file);
  save(save_file, ...
    'eigvals_sorted', 'cum_energy', 'eigvecs', 'sort_idx', ...
    'mean_vec', 'n99', 'rows_liq', 'cols_liq', 'frames', ...
    'case_name', '-v7.3');
  
  fprintf('[%s]  Done (%.1f s elapsed).\n\n', case_name, toc);
end

%% ── Slide 5 figure ───────────────────────────────────────────────────
fprintf('\nPlotting...\n');

fig = figure('Color', 'white', 'Position', [100 100 1100 620]);
ax  = axes(fig);
hold on;

n_cases = size(CASES, 1);
mode_90_vals = zeros(n_cases, 1);

for c = 1:n_cases
  cum = cum_energy_all{c};
  N   = numel(cum);
  
  % Log-log: skip zero values for log scale
  modes = (1:N)';
  plot(modes, cum * 100, ...
    'Color',     COLOURS(c,:), ...
    'LineWidth', LINE_W(c), ...
    'DisplayName', CASES{c,2});
  
  % Record where this case hits 90%
  idx90 = find(cum >= 0.90, 1);
  if ~isempty(idx90)
    mode_90_vals(c) = idx90;
    plot(idx90, 90, 'o', ...
      'Color',           COLOURS(c,:), ...
      'MarkerFaceColor', COLOURS(c,:), ...
      'MarkerSize',      8, ...
      'HandleVisibility','off');
  end
end

%% ── Axes scaling & labels ────────────────────────────────────────────
N_max = max(cellfun(@numel, cum_energy_all));
set(ax, ...
  'XScale',              'log', ...
  'YScale',              'linear', ...
  'XLim',                [1, N_max], ...
  'YLim',                [0, 100], ...
  'FontSize',            18, ...
  'FontName',            'Helvetica', ...
  'TickDir',             'in', ...
  'TickLabelInterpreter','latex', ...
  'LineWidth',           1.4, ...
  'Box',                 'on', ...
  'XMinorTick',          'on', ...
  'XColor',              TEXT, ...
  'YColor',              TEXT);

xlabel('Mode number (log scale)', ...
  'Interpreter', 'latex', 'FontSize', 24);
ylabel('Cumulative energy (\%)', ...
  'Interpreter', 'latex', 'FontSize', 24);
title('POD energy convergence liquid phase', ...
  'Interpreter', 'latex', 'FontSize', 20);

grid on;
ax.GridAlpha      = 0.15;
ax.GridColor      = MUTED;
ax.MinorGridAlpha = 0.08;

%% ── Threshold lines ──────────────────────────────────────────────────
THRESH_LABELS = {'80\%', '90\%'};
THRESH_VALS   = [80, 90];
for t = 1:numel(THRESH_VALS)
  yline(THRESH_VALS(t), '--', ...
    'Color',     MUTED, ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
  text(1.3, THRESH_VALS(t) + 1.5, ...
    [THRESH_LABELS{t} ' energy'], ...
    'Interpreter', 'latex', ...
    'FontSize',    16, ...
    'Color',       MUTED, ...
    'VerticalAlignment', 'bottom');
end

%% ── Key number callout (582 modes / 90%) ─────────────────────────────
xline(582, ':', ...
  'Color',     TEXT, ...
  'LineWidth', 1.0, ...
  'HandleVisibility', 'off');
text(582 * 1.05, 30, '582 modes', ...
  'Interpreter', 'latex', ...
  'FontSize',    15, ...
  'Color',       TEXT, ...
  'FontWeight',  'bold');

%% ── Noise region shading ─────────────────────────────────────────────
patch([1500 N_max N_max 1500], [0 0 100 100], ...
  [0.95 0.85 0.85], ...
  'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
  'HandleVisibility', 'off');
text(1500 * 1.05, 50, 'Noise', ...
  'Interpreter',       'latex', ...
  'FontSize',           16, ...
  'Color',              NOISE_COL(), ...
  'VerticalAlignment', 'middle');

%% ── G3 highlight annotation ─────────────────────────────────────────
cum_g3  = cum_energy_all{G3_IDX};
idx90g3 = find(cum_g3 >= 0.90, 1);
if ~isempty(idx90g3)
  annotation(fig, 'textarrow', ...
    ax_norm(ax, idx90g3 * 2.0, 80), ...
    ax_norm(ax, idx90g3,       90), ...
    'String',      'G3: slowest convergence', ...
    'Interpreter', 'latex', ...
    'FontSize',    15, ...
    'Color',       COLOURS(G3_IDX,:), ...
    'HeadStyle',   'vback2', ...
    'HeadLength',  8, ...
    'HeadWidth',   6, ...
    'LineWidth',   1.2);
end

%% ── Legend ───────────────────────────────────────────────────────────
leg = legend('Location', 'southeast', ...
  'Interpreter', 'latex', ...
  'FontSize',    18, ...
  'Box',         'on');
leg.EdgeColor = MUTED;

%% ── Save ─────────────────────────────────────────────────────────────
out_path = fullfile(OUT_DIR, 'energy_spectrum.png');
exportgraphics(fig, out_path, 'Resolution', 300, 'BackgroundColor', 'white');
fprintf('Saved:  %s\n', out_path);
fprintf('Total elapsed: %.1f s\n', toc);

%% ── Local helpers ────────────────────────────────────────────────────
function c = NOISE_COL()
c = [0.7 0.1 0.1];
end

% Convert data coordinates → normalised figure coords for annotation
function n = ax_norm(ax, x, y)
xl = get(ax, 'XLim'); yl = get(ax, 'YLim');
pos = ax.Position;
% x is on log scale
if strcmp(ax.XScale, 'log')
  nx = (log10(x) - log10(xl(1))) / (log10(xl(2)) - log10(xl(1)));
else
  nx = (x - xl(1)) / (xl(2) - xl(1));
end
ny = (y - yl(1)) / (yl(2) - yl(1));
n  = [pos(1) + nx * pos(3), pos(2) + ny * pos(4)];
end
