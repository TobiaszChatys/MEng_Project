function captureFrames(S, frames, filename)

[~, name, ~] = fileparts(filename);
aviName = [name, '_animation.avi'];

% 16:9 resolution: 1920x1080 (Full HD)
fig = figure("Position", [100 100 1920 1080], "color", "w", "Visible", "off");

try
  v = VideoWriter(aviName, 'Motion JPEG AVI');
  v.FrameRate = 8;
  v.Quality = 95;
  open(v);

  fprintf('Starting animation capture: %s\n', aviName);
  animateVelocityVectors(S, frames, v, fig);

  close(v);
  fprintf('AVI saved to: %s\n', aviName);

catch ME
  close(fig);
  fprintf('Error during video capture: %s\n', ME.message);
  rethrow(ME);
end

close(fig);
end
