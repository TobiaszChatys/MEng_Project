function animateVelocityVectors(S, frames, videoWriter, fig)

if nargin < 4 || isempty(fig)
  fig = gcf;
end

isCapturing = exist('videoWriter', 'var') && ~isempty(videoWriter);

% ------------------------------------------------------------------ %
%  Pre-compute global color limits from the raw struct arrays so the  %
%  colorbar scale stays fixed across all frames.                      %
% ------------------------------------------------------------------ %
n = min(frames, size(S.all_u_matrix_air, 3));
fprintf('Computing global velocity ranges across %d frames...\n', n);

% Pull raw velocity arrays directly -- no masking, no modification
Z_air_all = hypot(S.all_u_matrix_air(:,:,1:n), S.all_v_matrix_air(:,:,1:n));
Z_liq_all = hypot(S.all_u_matrix_liquid(:,:,1:n), S.all_v_matrix_liquid(:,:,1:n));

clim_air    = [0, prctile(Z_air_all(:), 99)];
clim_liquid = [0, prctile(Z_liq_all(:), 99)];
fprintf('CLim air: [%.3f, %.3f] m/s | liquid: [%.3f, %.3f] m/s\n', ...
  clim_air(1), clim_air(2), clim_liquid(1), clim_liquid(2));

clear Z_air_all Z_liq_all;

% Initialize graphics handles
ax1 = []; ax2 = []; h = struct();

for frame = 1:n
  if ~isvalid(fig)
    break;
  end
  
  % ---- Raw data: direct struct access, unit conversion only ---- %
  % Film height (PLIF) -- remove unphysical negative heights so
  % createMasks can interpolate cleanly; no velocity data is changed.
  X3 = S.film_x_position(1, :) * 1e3;
  Y3 = S.smoothed_film_height_matrix_out(:, frame) * 1e3;
  valid = Y3 >= 0;
  X3 = X3(valid);  Y3 = Y3(valid);
  
  % Air PIV -- raw velocities, NaNs zeroed for clean presentation
  X1 = S.all_transposed_x_position_matrix_air(:,:,frame) * 1e3;
  Y1 = S.all_transposed_y_position_matrix_air(:,:,frame) * 1e3;
  U1 = S.all_u_matrix_air(:,:,frame);
  V1 = S.all_v_matrix_air(:,:,frame);
  nan_air = isnan(U1) | isnan(V1);
  U1(nan_air) = 0; V1(nan_air) = 0;
  Z1 = hypot(U1, V1);
  
  % Liquid PIV -- raw velocities, NaNs zeroed for clean presentation
  X2 = S.all_transposed_x_position_matrix_liquid(:,:,frame) * 1e3;
  Y2 = S.all_transposed_y_position_matrix_liquid(:,:,frame) * 1e3;
  U2 = S.all_u_matrix_liquid(:,:,frame);
  V2 = S.all_v_matrix_liquid(:,:,frame);
  nan_liq = isnan(U2) | isnan(V2);
  U2(nan_liq) = 0; V2(nan_liq) = 0;
  Z2 = hypot(U2, V2);
  % -------------------------------------------------------------- %
  
  [Z1_masked, Z2_masked] = createMasks(X1, Y1, X2, Y2, Z1, Z2, X3, Y3);
  
  [ax1, ax2, h] = plotVectorMap(X1, Y1, U1, V1, X2, Y2, U2, V2, X3, Y3, ...
    Z1_masked, Z2_masked, ax1, ax2, h, clim_air, clim_liquid);
  
  title(ax2, sprintf('Velocity Vector Map: Frame %d / %d', frame, n), ...
    'FontSize', 32, 'FontWeight', 'bold', 'Interpreter', 'latex');
  
  if isCapturing
    drawnow;  % flush rendering queue so the frame is fully painted before capture
    frameCapture = getframe(fig);
    writeVideo(videoWriter, frameCapture);
  else
    drawnow limitrate;
    pause(0.05);
  end
end
end
