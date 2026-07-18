%% slide6_mode_reconstruction.m
%  Produces Slide 6 — the "Killer Animation": Raw PIV (static, top)
%  vs POD reconstruction with increasing mode counts (animated, bottom).
%  Presentation-ready 16:9 layout for embedding in PowerPoint.
%
%  Run from MEng_Project root:
%    >> run scripts/slide6_mode_reconstruction.m

tic; clc; close all;
addpath(genpath('src'));

%% ── Settings ─────────────────────────────────────────────────────────
CASE_NAME   = 'L8_G9';
FRAME_IDX   = 107;
MODE_COUNTS = [5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 500, 800, 953];
HOLD_SEC    = 1.5;                 % seconds per step
FPS         = 10;
FIG_SIZE    = [100 100 1920 1080]; % 16:9 Full HD

SCRIPT_DIR  = fileparts(mfilename('fullpath'));
PROJ_ROOT   = fileparts(SCRIPT_DIR);
DATA_DIR    = fullfile(PROJ_ROOT, 'data', 'Cases');
POD_DIR     = fullfile(SCRIPT_DIR, 'Results', 'POD_data');
OUT_DIR     = fullfile(SCRIPT_DIR, 'slide6_output');
if ~exist(POD_DIR, 'dir'), mkdir(POD_DIR); end
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 1 — Load raw data
%% ═══════════════════════════════════════════════════════════════════════
fprintf('Loading raw data: %s.mat ...\n', CASE_NAME);
mat_file = fullfile(DATA_DIR, [CASE_NAME '.mat']);
tmp = load(mat_file);
S   = tmp.S2P_PIV_full_mat_vars;
clear tmp

[rows_liq, cols_liq, N_frames] = size(S.all_u_matrix_liquid);
sp = rows_liq * cols_liq;

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 2 — Load POD results
%% ═══════════════════════════════════════════════════════════════════════
pod_file = fullfile(POD_DIR, ['POD_Results_' CASE_NAME '.mat']);
if ~exist(pod_file, 'file')
  error('POD results not found at %s. Please run POD analysis first.', pod_file);
end

fprintf('Loading POD results...\n');
d = load(pod_file);

if isfield(d, 'Spatial_modes'), spatial_modes = double(d.Spatial_modes);
elseif isfield(d, 'spatial_modes'), spatial_modes = double(d.spatial_modes);
else, spatial_modes = []; end

if isfield(d, 'mean_matrix'), mean_vec = double(d.mean_matrix(:));
elseif isfield(d, 'mean_vec'), mean_vec = double(d.mean_vec(:));
else, mean_vec = []; end

if isfield(d, 'cumulative_energy'), cum_energy = double(d.cumulative_energy(:));
elseif isfield(d, 'cum_energy'), cum_energy = double(d.cum_energy(:));
else, error('Cumulative energy not found in POD file'); end

% Load eigenvectors if spatial modes are missing
eigvecs = []; sort_idx = [];
if isempty(spatial_modes)
  if isfield(d, 'eigenvectors_matrix'), eigvecs = double(d.eigenvectors_matrix);
  elseif isfield(d, 'eigvecs'), eigvecs = double(d.eigvecs); end
  if isfield(d, 'sort_index'), sort_idx = double(d.sort_index(:));
  elseif isfield(d, 'sort_idx'), sort_idx = double(d.sort_idx(:)); end
end

clear d

% ── Compute spatial modes if not cached ───────────────────────────────
MAX_MODES_NEEDED = max(MODE_COUNTS);
if isempty(spatial_modes)
  fprintf('  Spatial modes not cached — computing from eigenvectors...\n');
  if isempty(eigvecs) || isempty(sort_idx)
    error('Neither Spatial_modes nor eigenvectors found in %s', pod_file);
  end
  
  fprintf('  Building snapshot matrix (%d frames)...\n', N_frames);
  U_vec = zeros(sp, N_frames);
  V_vec = zeros(sp, N_frames);
  for bfr = 1:N_frames
    U_vec(:,bfr) = reshape(S.all_u_matrix_liquid(:,:,bfr), [], 1);
    V_vec(:,bfr) = reshape(S.all_v_matrix_liquid(:,:,bfr), [], 1);
  end
  Snap = [U_vec; V_vec];
  Snap(isnan(Snap)) = 0;
  clear U_vec V_vec
  
  if isempty(mean_vec)
    mean_vec = mean(Snap, 2);
  end
  Snap_fluc = Snap - mean_vec;
  clear Snap
  
  K_COMPUTE = min(MAX_MODES_NEEDED, size(eigvecs, 2));
  fprintf('  Computing %d spatial modes...\n', K_COMPUTE);
  V_sorted = eigvecs(:, sort_idx(1:K_COMPUTE));
  spatial_modes = Snap_fluc * V_sorted;
  for k = 1:K_COMPUTE
    spatial_modes(:,k) = spatial_modes(:,k) / norm(spatial_modes(:,k));
  end
  clear Snap_fluc V_sorted eigvecs sort_idx
  
  % Cache for next time
  fprintf('  Appending Spatial_modes to %s ...\n', pod_file);
  Spatial_modes = spatial_modes;                     %#ok
  mean_matrix   = mean_vec;                          %#ok
  save(pod_file, 'Spatial_modes', 'mean_matrix', '-append');
  clear Spatial_modes mean_matrix
  fprintf('  Spatial modes cached.\n');
end

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 3 — Prepare Data
%% ═══════════════════════════════════════════════════════════════════════
% Coordinates
X_SHIFT = 11;
X = S.all_transposed_x_position_matrix_liquid(:,:,FRAME_IDX) * 1e3 + X_SHIFT;
Y = S.all_transposed_y_position_matrix_liquid(:,:,FRAME_IDX) * 1e3;

% Raw Velocity
U_raw = S.all_u_matrix_liquid(:,:,FRAME_IDX);
V_raw = S.all_v_matrix_liquid(:,:,FRAME_IDX);

% Interface
X_film = S.film_x_position(1,:) * 1e3 + X_SHIFT;
Y_film = S.smoothed_film_height_matrix_out(:, FRAME_IDX)' * 1e3;
valid_film = Y_film >= 0;
X_film = X_film(valid_film);
Y_film = Y_film(valid_film);

% POD Reconstruction Setup
nan_mask  = isnan(U_raw) | isnan(V_raw);
U_raw(nan_mask) = 0;
V_raw(nan_mask) = 0;
Z_raw = hypot(U_raw, V_raw);
snap_raw  = [reshape(U_raw, [], 1); reshape(V_raw, [], 1)];
snap_raw(isnan(snap_raw)) = 0;
snap_fluc = snap_raw - mean_vec;

% Interface mask — zero out anything at or above the film height
Y_film_interp = interp1(X_film, Y_film, X(1,:), 'linear', 'extrap');
mask_above    = Y >= repmat(Y_film_interp, size(Y,1), 1);

Z_raw(mask_above)  = 0;
U_raw(mask_above)  = 0;
V_raw(mask_above)  = 0;

% Calculate max modes available
max_modes_avail = size(spatial_modes, 2);
MODE_COUNTS = MODE_COUNTS(MODE_COUNTS <= max_modes_avail);

% Pre-calculate coefficients
a_coeffs = spatial_modes' * snap_fluc;

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 4 — Render Video
%% ═══════════════════════════════════════════════════════════════════════
fprintf('Rendering video...\n');

avi_path = fullfile(OUT_DIR, [CASE_NAME '_mode_reconstruction.avi']);
mp4_path = fullfile(OUT_DIR, [CASE_NAME '_mode_reconstruction.mp4']);

v = VideoWriter(avi_path, 'Motion JPEG AVI');
v.FrameRate = FPS;
v.Quality   = 95;
open(v);

fig = figure('Position', FIG_SIZE, 'Color', 'w', 'Visible', 'off');

% Visual Settings
CLIM_VAL   = [0, prctile(Z_raw(:), 99)];
FILM_COL   = [0, 1, 1];
FILM_LW    = 4;
FONT_TITLE = 28;
FONT_PANEL = 24;
FONT_LABEL = 22;
FONT_AXIS  = 18;
YLIM_MAX   = 4;           % mm — focus on liquid film region

% ── Layout: 2 stacked subplots with manual positioning ────────────────
left_m   = 0.07;
right_m  = 0.10;           % space for individual colorbars
top_m    = 0.06;
bottom_m = 0.08;
gap_v    = 0.05;
cb_w     = 0.02;
cb_gap   = 0.01;

avail_w  = 1 - left_m - right_m;
avail_h  = 1 - top_m - bottom_m - gap_v;
panel_h  = avail_h / 2;

% ── Top Plot: Static Raw PIV ──────────────────────────────────────────
ax1 = axes(fig, 'Position', [left_m, bottom_m + panel_h + gap_v, avail_w, panel_h]);
hold(ax1, 'on');
pcolor(ax1, X, Y, Z_raw); shading(ax1, 'flat');
quiver(ax1, X, Y, U_raw, V_raw, 1.2, 'w');
plot(ax1, X_film, Y_film, 'Color', FILM_COL, 'LineWidth', FILM_LW);

colormap(ax1, winter);
clim(ax1, CLIM_VAL);
xlim(ax1, [-15, 15]);
ylim(ax1, [0, YLIM_MAX]);
set(ax1, ...
  'XTick',               -15:5:15, ...
  'XTickLabel',          [], ...
  'FontSize',            FONT_AXIS, ...
  'TickDir',             'out', ...
  'TickLabelInterpreter','latex', ...
  'LineWidth',           1.5, ...
  'Layer',               'top');
ylabel(ax1, '$Y$ (mm)', 'FontSize', FONT_LABEL, 'Interpreter', 'latex');

% Colorbar for top panel
cb1_left = left_m + avail_w + cb_gap;
cb1_bot  = bottom_m + panel_h + gap_v;
cb1 = colorbar(ax1, 'Position', [cb1_left, cb1_bot, cb_w, panel_h]);
cb1.Label.String         = '$U$ (m s$^{-1}$)';
cb1.Label.Interpreter    = 'latex';
cb1.Label.FontSize       = FONT_AXIS;
cb1.FontSize             = FONT_AXIS - 2;
cb1.TickLabelInterpreter = 'latex';

% Panel label
xl = xlim(ax1);  yl = ylim(ax1);
text(ax1, xl(1) + 0.5, yl(2) - 0.30, 'Raw PIV', ...
  'FontSize',          FONT_PANEL, ...
  'FontWeight',        'bold', ...
  'Color',             'w', ...
  'BackgroundColor',   [0.1 0.1 0.1], ...
  'VerticalAlignment', 'top', ...
  'Interpreter',       'none', ...
  'Margin',            6);

% ── Bottom Plot: Animated POD Reconstruction ──────────────────────────
ax2 = axes(fig, 'Position', [left_m, bottom_m, avail_w, panel_h]);
hold(ax2, 'on');

hm = pcolor(ax2, X, Y, zeros(size(Z_raw))); shading(ax2, 'flat');
hq = quiver(ax2, X, Y, zeros(size(U_raw)), zeros(size(V_raw)), 1.2, 'w');
hf = plot(ax2, X_film, Y_film, 'Color', FILM_COL, 'LineWidth', FILM_LW);

colormap(ax2, winter);
clim(ax2, CLIM_VAL);
xlim(ax2, [-15, 15]);
ylim(ax2, [0, YLIM_MAX]);
set(ax2, ...
  'XTick',               -15:5:15, ...
  'FontSize',            FONT_AXIS, ...
  'TickDir',             'out', ...
  'TickLabelInterpreter','latex', ...
  'LineWidth',           1.5, ...
  'Layer',               'top');
xlabel(ax2, '$X$ (mm)', 'FontSize', FONT_LABEL, 'Interpreter', 'latex');
ylabel(ax2, '$Y$ (mm)', 'FontSize', FONT_LABEL, 'Interpreter', 'latex');

% Colorbar for bottom panel
cb2_left = left_m + avail_w + cb_gap;
cb2 = colorbar(ax2, 'Position', [cb2_left, bottom_m, cb_w, panel_h]);
cb2.Label.String         = '$U$ (m s$^{-1}$)';
cb2.Label.Interpreter    = 'latex';
cb2.Label.FontSize       = FONT_AXIS;
cb2.FontSize             = FONT_AXIS - 2;
cb2.TickLabelInterpreter = 'latex';

% Mode-count overlay label (updated each step)
xl2 = xlim(ax2);  yl2 = ylim(ax2);
h_lbl = text(ax2, xl2(1) + 0.5, yl2(2) - 0.30, '', ...
  'FontSize',          FONT_PANEL, ...
  'FontWeight',        'bold', ...
  'Color',             'w', ...
  'BackgroundColor',   [0.1 0.1 0.1], ...
  'VerticalAlignment', 'top', ...
  'Interpreter',       'none', ...
  'Margin',            6);

linkaxes([ax1, ax2], 'xy');

% ── Super-title ───────────────────────────────────────────────────────
h_stitle = annotation(fig, 'textbox', ...
  [0.05, 1 - top_m + 0.005, 0.9, top_m - 0.01], ...
  'String', sprintf('POD Mode Reconstruction %s Frame %d', ...
  CASE_NAME, FRAME_IDX), ...
  'FontSize',            FONT_TITLE, ...
  'FontWeight',          'bold', ...
  'HorizontalAlignment', 'center', ...
  'VerticalAlignment',   'middle', ...
  'EdgeColor',           'none', ...
  'Interpreter',         'none');

% ── Animation Loop ────────────────────────────────────────────────────
frames_per_step = round(HOLD_SEC * FPS);

for k = MODE_COUNTS
  % Reconstruct
  snap_recon = mean_vec + spatial_modes(:,1:k) * a_coeffs(1:k);
  U_rec = reshape(snap_recon(1:sp), rows_liq, cols_liq);
  V_rec = reshape(snap_recon(sp+1:end), rows_liq, cols_liq);
  
  % Mask — original NaNs + zero out above interface
  U_rec(nan_mask)    = 0;
  V_rec(nan_mask)    = 0;
  U_rec(mask_above)  = 0;
  V_rec(mask_above)  = 0;
  Z_rec = hypot(U_rec, V_rec);
  
  % Update Plot
  set(hm, 'CData', Z_rec);
  set(hq, 'UData', U_rec, 'VData', V_rec);
  
  % Update overlay label
  energy_pct = cum_energy(k) * 100;
  label_str = sprintf('%d Modes - %.1f%% Energy', k, energy_pct);
  set(h_lbl, 'String', label_str);
  uistack(h_lbl, 'top');
  
  % Capture
  drawnow;
  frame = getframe(fig);
  for i = 1:frames_per_step
    writeVideo(v, frame);
  end
  fprintf('  Processed: %s\n', label_str);
end

close(v);
fprintf('AVI written: %s\n', avi_path);

%% ── Save final frame as .fig for later editing ──────────────────────
fig_path = fullfile(OUT_DIR, [CASE_NAME '_mode_reconstruction_final.fig']);
set(fig, 'Visible', 'on');
savefig(fig, fig_path);
fprintf('Figure saved: %s\n', fig_path);

%% ── Save final frame as high-res PNG ─────────────────────────────────
png_path = fullfile(OUT_DIR, [CASE_NAME '_mode_reconstruction_final.png']);
exportgraphics(fig, png_path, 'Resolution', 300, 'BackgroundColor', 'white');
fprintf('PNG saved: %s\n', png_path);

close(fig);

%% ── ffmpeg: AVI → MP4 ───────────────────────────────────────────────
fprintf('Converting to MP4 via ffmpeg...\n');
cmd = sprintf(['env -u LD_LIBRARY_PATH ffmpeg -y -i "%s" ' ...
  '-c:v libx264 -preset slow -crf 18 -pix_fmt yuv420p ' ...
  '"%s" 2>&1'], avi_path, mp4_path);
[status, result] = system(cmd);
if status == 0
  fprintf('MP4 saved: %s\n', mp4_path);
  delete(avi_path);
else
  fprintf('ffmpeg failed:\n%s\n', result);
end

fprintf('\nDone! Total elapsed: %.1f s\n', toc);
