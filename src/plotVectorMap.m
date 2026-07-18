function [ax1, ax2, h] = plotVectorMap(X1, Y1, U1, V1, X2, Y2, U2, V2, X3, Y3, Z1_masked, Z2_masked, ax1, ax2, h, clim_air, clim_liquid)

% --- Configuration ---
AXIS_FONT_SIZE  = 24;
LABEL_FONT_SIZE = 30;
INSET_FONT_SIZE = 16;

% Coordinate shift to align PIV and PLIF grids
X_SHIFT = 11;

% Cyan [0 1 1]: absent from hot (no blue channel) and from winter (winter maxes at [0,1,0.5] -- B never reaches 1 outside its minimum)
INTERFACE_COLOR = [0, 1, 1];

% Zoom-inset region (data coordinates, mm)
ZOOM_XLIM = [-5, 5];
ZOOM_YLIM = [0,  4];

% Default CLim if not provided
if nargin < 16 || isempty(clim_air),    clim_air    = [0 1]; end
if nargin < 17 || isempty(clim_liquid), clim_liquid = [0 1]; end

% Ensure X3/Y3 are row vectors
X3_row = (X3(:))';
Y3_row = (Y3(:))';

% --- Layout (normalised figure units) ---
AX_POS    = [0.08, 0.12, 0.74, 0.80];
CB1_POS   = [0.86, 0.54, 0.025, 0.38];   % air   (hot)
CB2_POS   = [0.86, 0.12, 0.025, 0.38];   % liquid (winter)
INSET_POS = [0.10, 0.50, 0.25, 0.25];    % zoom inset -- 16:9 on a 16:9 figure (equal norm units)

% ------------------------------------------------------------------ %
%  Air Axis (background -- hot heatmap)                               %
% ------------------------------------------------------------------ %
if nargin < 13 || isempty(ax1) || ~isvalid(ax1)
  ax1 = axes('Position', AX_POS);
  hold(ax1, 'on');
  h.p1 = pcolor(ax1, X1 + X_SHIFT, Y1, Z1_masked);
  shading(ax1, 'flat');
  colormap(ax1, hot);
  clim(ax1, clim_air);
  set(ax1, 'CLimMode', 'manual');
  ax1.XAxis.Visible = 'off';
  ax1.YAxis.Visible = 'off';
  
  h.c1 = colorbar(ax1, 'Position', CB1_POS);
  h.c1.Label.String      = 'Air (m s$^{-1}$)';
  h.c1.Label.Interpreter = 'latex';
  h.c1.Label.FontSize    = AXIS_FONT_SIZE;
  h.c1.FontSize          = AXIS_FONT_SIZE - 4;
  h.c1.TickLabelInterpreter = 'latex';
else
  set(h.p1, 'CData', Z1_masked);
end

% ------------------------------------------------------------------ %
%  Liquid Axis (overlay -- winter heatmap + quivers + film height)    %
% ------------------------------------------------------------------ %
if nargin < 14 || isempty(ax2) || ~isvalid(ax2)
  ax2 = axes('Position', AX_POS, 'Color', 'none');
  hold(ax2, 'on');
  
  h.p2 = pcolor(ax2, X2 + X_SHIFT, Y2, Z2_masked);
  shading(ax2, 'flat');
  colormap(ax2, winter);
  clim(ax2, clim_liquid);
  set(ax2, 'CLimMode', 'manual');
  
  h.c2 = colorbar(ax2, 'Position', CB2_POS);
  h.c2.Label.String      = 'Liquid (m s$^{-1}$)';
  h.c2.Label.Interpreter = 'latex';
  h.c2.Label.FontSize    = AXIS_FONT_SIZE;
  h.c2.FontSize          = AXIS_FONT_SIZE - 4;
  h.c2.TickLabelInterpreter = 'latex';
  
  h.q1 = quiver(ax2, X1 + X_SHIFT, Y1, U1, V1, 0.8, 'k');
  h.q2 = quiver(ax2, X2 + X_SHIFT, Y2, U2, V2, 0.8, 'w');
  
  h.line = plot(ax2, X3_row + X_SHIFT, Y3_row, '-', ...
    'Color', INTERFACE_COLOR, 'LineWidth', 5);
  h.interface_label = text(ax2, mean(X3_row + X_SHIFT) + 4, mean(Y3_row) + 1.5, ...
    'Air--liquid interface', 'Color', INTERFACE_COLOR, 'FontSize', 24, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
  
  h.xline = xline(ax2, 0, 'k--', 'LineWidth', 2);
  
  linkaxes([ax1 ax2], 'xy');
  uistack(ax1, 'bottom');
  
  xlim(ax2, [-15, 15]);
  ylim(ax2, [0,  28]);
  set(ax2, 'XTick', -15:5:15, 'YTick', 0:2:28, ...
    'FontSize', AXIS_FONT_SIZE, ...
    'TickDir', 'out', ...
    'TickLabelInterpreter', 'latex', ...
    'LineWidth', 2.0, ...
    'XAxisLocation', 'bottom', 'YAxisLocation', 'left', ...
    'Layer', 'top');
  xlabel(ax2, '$X$ (mm)', 'FontSize', LABEL_FONT_SIZE, 'Interpreter', 'latex');
  ylabel(ax2, '$Y$ (mm)', 'FontSize', LABEL_FONT_SIZE, 'Interpreter', 'latex');
  
  % Dashed box on the main axes marking the zoomed region
  h.zoom_box = rectangle(ax2, ...
    'Position', [ZOOM_XLIM(1), ZOOM_YLIM(1), diff(ZOOM_XLIM), diff(ZOOM_YLIM)], ...
    'EdgeColor', INTERFACE_COLOR, 'LineStyle', '--', 'LineWidth', 2);
  
  % ------------------------------------------------------------ %
  %  Zoom inset -- liquid phase detail                            %
  % ------------------------------------------------------------ %
  h.ax_zoom = axes('Position', INSET_POS);
  hold(h.ax_zoom, 'on');
  
  h.zoom_p = pcolor(h.ax_zoom, X2 + X_SHIFT, Y2, Z2_masked);
  shading(h.ax_zoom, 'flat');
  colormap(h.ax_zoom, winter);
  clim(h.ax_zoom, clim_liquid);
  set(h.ax_zoom, 'CLimMode', 'manual');
  
  % Black arrows are visible against both winter colors and white background
  h.zoom_q2 = quiver(h.ax_zoom, X2 + X_SHIFT, Y2, U2, V2, 0.8, 'k');
  h.zoom_q1 = quiver(h.ax_zoom, X1 + X_SHIFT, Y1, U1, V1, 0.8, 'k', 'LineWidth', 0.5);
  
  h.zoom_line = plot(h.ax_zoom, X3_row + X_SHIFT, Y3_row, '-', ...
    'Color', INTERFACE_COLOR, 'LineWidth', 4);
  
  xlim(h.ax_zoom, ZOOM_XLIM);
  ylim(h.ax_zoom, ZOOM_YLIM);
  
  set(h.ax_zoom, ...
    'XTick', [], 'YTick', [], ...
    'Box', 'on','LineWidth', 2.5, ...
    'XColor', 'k', 'YColor', 'k', ...
    'Layer', 'top');
  
else
  % Update all graphics objects for subsequent frames
  set(h.p2,     'CData', Z2_masked);
  set(h.q1,     'XData', X1 + X_SHIFT, 'YData', Y1, 'UData', U1, 'VData', V1);
  set(h.q2,     'XData', X2 + X_SHIFT, 'YData', Y2, 'UData', U2, 'VData', V2);
  set(h.line,   'XData', X3_row + X_SHIFT, 'YData', Y3_row);
  set(h.interface_label, 'Position', [mean(X3_row + X_SHIFT) - 2, mean(Y3_row) + 1.5, 0]);
  uistack(h.interface_label, 'top');
  
  % Update inset
  set(h.zoom_p,  'CData', Z2_masked);
  set(h.zoom_q2, 'XData', X2 + X_SHIFT, 'YData', Y2, 'UData', U2, 'VData', V2);
  set(h.zoom_q1, 'XData', X1 + X_SHIFT, 'YData', Y1, 'UData', U1, 'VData', V1);
  set(h.zoom_line, 'XData', X3_row + X_SHIFT, 'YData', Y3_row);
end
end
