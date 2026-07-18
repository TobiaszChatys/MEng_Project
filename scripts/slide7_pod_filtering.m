%% slide7_pod_filtering.m
%  Produces Slide 7 — Before & After: POD Filtering Effect (Full Animation)
%
%  Two-panel animation (vertically stacked subplots):
%    TOP:    Raw PIV liquid phase — animated through all frames
%    BOTTOM: POD reconstruction using K_MODES most energetic modes
%
%  Both panels animate in sync.  Same layout/style as slide6.
%  Change K_MODES below to control filtering aggressiveness.
%
%  Run from MEng_Project root:
%    >> run scripts/slide7_pod_filtering.m
%
%  Output:  scripts/slide7_output/<CASE>_pod_filtering.mp4
%           scripts/slide7_output/<CASE>_pod_filtering_final.fig
%           scripts/slide7_output/<CASE>_pod_filtering_final.png

tic; clc; close all;
addpath(genpath('src'));

%% ── Settings ─────────────────────────────────────────────────────────
CASE_NAME  = 'L8_G6';             % 3D regime — clear wave structure
K_MODES    = 150;                 % ← number of most energetic modes to use
N_ANIM     = 300;                 % frames to animate
FPS        = 8;
FIG_SIZE   = [100 100 1920 1080]; % 16:9 Full HD
BLOCK_SIZE = 350;                 % frames per block (matches Incremental_POD)

SCRIPT_DIR = fileparts(mfilename('fullpath'));
PROJ_ROOT  = fileparts(SCRIPT_DIR);
DATA_DIR   = fullfile(PROJ_ROOT, 'data', 'Cases');
POD_DIR    = fullfile(SCRIPT_DIR, 'Results', 'POD_data');
OUT_DIR    = fullfile(SCRIPT_DIR, 'slide7_output');
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
N_ANIM = min(N_ANIM, N_frames);
fprintf('  Grid: %d x %d, %d frames (animating %d)\n', ...
  rows_liq, cols_liq, N_frames, N_ANIM);

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 2 — Load POD results
%% ═══════════════════════════════════════════════════════════════════════
pod_file = fullfile(POD_DIR, ['POD_Results_' CASE_NAME '.mat']);
if ~exist(pod_file, 'file')
  error('POD results not found at %s. Run slide5 or Incremental_POD first.', pod_file);
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
if isempty(spatial_modes)
  fprintf('  Spatial modes not cached — computing from eigenvectors...\n');
  if isempty(eigvecs) || isempty(sort_idx)
    error('Neither Spatial_modes nor eigenvectors found in %s', pod_file);
  end
  
  % Build fluctuation snapshot matrix
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
  
  K_COMPUTE = min(K_MODES, size(eigvecs, 2));
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

% Validate mode count
max_modes_avail = size(spatial_modes, 2);
if K_MODES > max_modes_avail
  warning('K_MODES (%d) > available modes (%d). Clamping.', K_MODES, max_modes_avail);
  K_MODES = max_modes_avail;
end
spatial_modes = spatial_modes(:, 1:K_MODES);
energy_pct = cum_energy(K_MODES) * 100;
fprintf('  Using %d modes (%.1f%% energy)\n', K_MODES, energy_pct);

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 3 — Pre-compute temporal coefficients for all animated frames
%% ═══════════════════════════════════════════════════════════════════════
fprintf('Pre-computing temporal coefficients for %d frames...\n', N_ANIM);
A_coeffs = zeros(K_MODES, N_ANIM);
for fr = 1:N_ANIM
  U2 = S.all_u_matrix_liquid(:,:,fr);
  V2 = S.all_v_matrix_liquid(:,:,fr);
  snap = [reshape(U2, [], 1); reshape(V2, [], 1)];
  snap(isnan(snap)) = 0;
  A_coeffs(:,fr) = spatial_modes' * (snap - mean_vec);
end
fprintf('  Coefficients computed: %.1f s\n', toc);

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 4 — Global colour limits (consistent across both panels/frames)
%% ═══════════════════════════════════════════════════════════════════════
fprintf('Computing colour limits...\n');
Z_liq_all = hypot(S.all_u_matrix_liquid(:,:,1:N_ANIM), ...
  S.all_v_matrix_liquid(:,:,1:N_ANIM));
CLIM_VAL = [0, prctile(Z_liq_all(:), 99)];
clear Z_liq_all
fprintf('  CLim liquid: [%.3f, %.3f] m/s\n', CLIM_VAL);

%% ═══════════════════════════════════════════════════════════════════════
%  PHASE 5 — Render animation (2 vertically-stacked subplots)
%% ═══════════════════════════════════════════════════════════════════════
fprintf('Rendering animation (%d frames @ %d FPS = %.1f s)...\n', ...
  N_ANIM, FPS, N_ANIM / FPS);

X_SHIFT    = 11;           % mm — PIV-PLIF grid alignment
FILM_COL   = [0, 1, 1];
FILM_LW    = 4;
FONT_TITLE = 28;
FONT_PANEL = 24;
FONT_LABEL = 22;
FONT_AXIS  = 18;
YLIM_MAX   = 6;            % mm — focus on liquid film region

panel_labels = { ...
  'Raw PIV', ...
  sprintf('%d Modes - %.1f%% Energy', K_MODES, energy_pct)};

% ── Layout: 2 stacked subplots with manual positioning ────────────────
left_m   = 0.07;
right_m  = 0.10;
top_m    = 0.06;
bottom_m = 0.08;
gap_v    = 0.05;
cb_w     = 0.02;
cb_gap   = 0.01;

avail_w  = 1 - left_m - right_m;
avail_h  = 1 - top_m - bottom_m - gap_v;
panel_h  = avail_h / 2;

avi_path = fullfile(OUT_DIR, [CASE_NAME '_pod_filtering.avi']);
mp4_path = fullfile(OUT_DIR, [CASE_NAME '_pod_filtering.mp4']);

vw = VideoWriter(avi_path, 'Motion JPEG AVI');
vw.FrameRate = FPS;
vw.Quality   = 95;
open(vw);

fig = figure('Position', FIG_SIZE, 'Color', 'w', 'Visible', 'off');

h = struct();   % graphics handles

for fr = 1:N_ANIM
  % ── Extract raw data for this frame ───────────────────────────────
  X = S.all_transposed_x_position_matrix_liquid(:,:,fr) * 1e3 + X_SHIFT;
  Y = S.all_transposed_y_position_matrix_liquid(:,:,fr) * 1e3;
  U_raw = S.all_u_matrix_liquid(:,:,fr);
  V_raw = S.all_v_matrix_liquid(:,:,fr);
  nan_mask = isnan(U_raw) | isnan(V_raw);
  U_raw(nan_mask) = 0;
  V_raw(nan_mask) = 0;
  
  % Film height (interface)
  X_film = S.film_x_position(1,:) * 1e3 + X_SHIFT;
  Y_film = S.smoothed_film_height_matrix_out(:,fr)' * 1e3;
  valid_film = Y_film >= 0;
  X_film = X_film(valid_film);
  Y_film = Y_film(valid_film);
  
  % Interface mask — zero out above film height
  Y_film_interp = interp1(X_film, Y_film, X(1,:), 'linear', 'extrap');
  mask_above    = Y >= repmat(Y_film_interp, size(Y,1), 1);
  
  Z_raw = hypot(U_raw, V_raw);
  Z_raw(mask_above)  = 0;
  U_raw(mask_above)  = 0;
  V_raw(mask_above)  = 0;
  
  % ── POD reconstruction ────────────────────────────────────────────
  a = A_coeffs(:,fr);
  snap_rec = mean_vec + spatial_modes * a;
  U_rec = reshape(snap_rec(1:sp),     rows_liq, cols_liq);
  V_rec = reshape(snap_rec(sp+1:end), rows_liq, cols_liq);
  U_rec(nan_mask)   = 0;
  V_rec(nan_mask)   = 0;
  U_rec(mask_above) = 0;
  V_rec(mask_above) = 0;
  Z_rec = hypot(U_rec, V_rec);
  
  if fr == 1
    % ── First frame: create all graphics objects ──────────────────
    
    % ── Top Plot: Raw PIV ─────────────────────────────────────────
    ax1 = axes(fig, 'Position', [left_m, bottom_m + panel_h + gap_v, avail_w, panel_h]);
    hold(ax1, 'on');
    h.hm1 = pcolor(ax1, X, Y, Z_raw); shading(ax1, 'flat');
    h.q1  = quiver(ax1, X, Y, U_raw, V_raw, 1.2, 'w');
    h.f1  = plot(ax1, X_film, Y_film, 'Color', FILM_COL, 'LineWidth', FILM_LW);
    
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
    h.lbl1 = text(ax1, xl(1) + 0.5, yl(2) - 0.30, panel_labels{1}, ...
      'FontSize',          FONT_PANEL, ...
      'FontWeight',        'bold', ...
      'Color',             'w', ...
      'BackgroundColor',   [0.1 0.1 0.1], ...
      'VerticalAlignment', 'top', ...
      'Interpreter',       'none', ...
      'Margin',            6);
    
    % ── Bottom Plot: POD Reconstruction ───────────────────────────
    ax2 = axes(fig, 'Position', [left_m, bottom_m, avail_w, panel_h]);
    hold(ax2, 'on');
    h.hm2 = pcolor(ax2, X, Y, Z_rec); shading(ax2, 'flat');
    h.q2  = quiver(ax2, X, Y, U_rec, V_rec, 1.2, 'w');
    h.f2  = plot(ax2, X_film, Y_film, 'Color', FILM_COL, 'LineWidth', FILM_LW);
    
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
    
    % Panel label
    xl2 = xlim(ax2);  yl2 = ylim(ax2);
    h.lbl2 = text(ax2, xl2(1) + 0.5, yl2(2) - 0.30, panel_labels{2}, ...
      'FontSize',          FONT_PANEL, ...
      'FontWeight',        'bold', ...
      'Color',             'w', ...
      'BackgroundColor',   [0.1 0.1 0.1], ...
      'VerticalAlignment', 'top', ...
      'Interpreter',       'none', ...
      'Margin',            6);
    
    linkaxes([ax1, ax2], 'xy');
    
    % ── Super-title ──────────────────────────────────────────────
    h.stitle = annotation(fig, 'textbox', ...
      [0.05, 1 - top_m + 0.005, 0.9, top_m - 0.01], ...
      'String', sprintf( ...
      'POD Filtering Effect: %s Frame %d / %d', ...
      CASE_NAME, fr, N_ANIM), ...
      'FontSize',            FONT_TITLE, ...
      'FontWeight',          'bold', ...
      'HorizontalAlignment', 'center', ...
      'VerticalAlignment',   'middle', ...
      'EdgeColor',           'none', ...
      'Interpreter',         'none');
  else
    % ── Subsequent frames: update existing objects ────────────────
    % Top panel (raw)
    set(h.hm1, 'CData', Z_raw);
    set(h.q1,  'XData', X, 'YData', Y, 'UData', U_raw, 'VData', V_raw);
    set(h.f1,  'XData', X_film, 'YData', Y_film);
    uistack(h.lbl1, 'top');
    
    % Bottom panel (reconstruction)
    set(h.hm2, 'CData', Z_rec);
    set(h.q2,  'XData', X, 'YData', Y, 'UData', U_rec, 'VData', V_rec);
    set(h.f2,  'XData', X_film, 'YData', Y_film);
    uistack(h.lbl2, 'top');
    
    set(h.stitle, 'String', ...
      sprintf('POD Filtering Effect %s Frame %d / %d', ...
      CASE_NAME, fr, N_ANIM));
  end
  
  drawnow;
  frame_capture = getframe(fig);
  writeVideo(vw, frame_capture);
  
  if mod(fr, 25) == 0
    fprintf('  Frame %d / %d  (%.1f s)\n', fr, N_ANIM, toc);
  end
end

close(vw);
fprintf('AVI written: %s\n', avi_path);

%% ── Save final frame as .fig for later editing ──────────────────────
fig_path = fullfile(OUT_DIR, [CASE_NAME '_pod_filtering_final.fig']);
set(fig, 'Visible', 'on');
savefig(fig, fig_path);
fprintf('Figure saved: %s\n', fig_path);

%% ── Save final frame as high-res PNG ─────────────────────────────────
png_path = fullfile(OUT_DIR, [CASE_NAME '_pod_filtering_final.png']);
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
  fprintf('AVI removed.\n');
else
  fprintf('ffmpeg conversion failed:\n%s\n', result);
  fprintf('AVI kept at: %s\n', avi_path);
end

fprintf('\nDone! Total elapsed: %.1f s\n', toc);
