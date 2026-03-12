function [ax1, ax2, h] = plotVectorMap(X1, Y1, U1, V1, X2, Y2, U2, V2, X3, Y3, Z1_masked, Z2_masked, ax1, ax2, h)

% Air Axis (Bottom)
if nargin < 13 || isempty(ax1) || ~isvalid(ax1)
  ax1 = axes("Position", [0.1, 0.1, 0.8, 0.8]);
  c1 = colorbar(ax1, "Position", [0.92, 0.1, 0.02, 0.8]);
  c1.Label.String = 'Air Velocity Magnitude (m/s)';
  hold(ax1, 'on');
  h.p1 = pcolor(ax1, X1 + 11, Y1, Z1_masked);
  shading(ax1, 'flat');
  colormap(ax1, hot);
else
  % Update only CData for pcolor to keep it 2D and below quivers
  set(h.p1, 'CData', Z1_masked);
end

% Liquid Axis (Top)
if nargin < 14 || isempty(ax2) || ~isvalid(ax2)
  ax2 = axes('Position', [0.1, 0.1, 0.8, 0.8], 'Color', 'none');
  c2 = colorbar(ax2, "Position", [0.96, 0.1, 0.02, 0.8]);
  c2.Label.String = 'Liquid Velocity Magnitude (m/s)';
  hold(ax2, 'on');
  
  % Liquid heatmap
  h.p2 = pcolor(ax2, X2 + 11, Y2, Z2_masked);
  shading(ax2, 'flat');
  colormap(ax2, jet);
  
  % Put BOTH quivers on the top axis to ensure they are always visible
  h.q1 = quiver(ax2, X1 + 11, Y1, U1, V1, 0.8, 'k');
  h.q2 = quiver(ax2, X2 + 11, Y2, U2, V2, 0.8, 'w');
  
  h.line = plot(ax2, X3 + 11, Y3, 'g', 'LineWidth', 4);
  h.xline = xline(ax2, 0, 'k--', 'LineWidth', 2);
  
  linkaxes([ax1 ax2], 'xy');
  uistack(ax1, 'bottom');
  
  xlim(ax2, [-10 10]);
  ylim(ax2, [0 10]);
  xlabel(ax2, 'X Position (mm)');
  ylabel(ax2, 'Y Position (mm)');
  y_ticks = 0:2:28;
  set(ax2, 'YTick', y_ticks, 'YTickLabel', y_ticks);
else
  % Update liquid heatmap
  set(h.p2, 'CData', Z2_masked);
  
  % Update BOTH quivers
  set(h.q1, 'XData', X1 + 11, 'YData', Y1, 'UData', U1, 'VData', V1);
  set(h.q2, 'XData', X2 + 11, 'YData', Y2, 'UData', U2, 'VData', V2);
  
  % Update interface line
  set(h.line, 'XData', X3 + 11, 'YData', Y3);
end
end
