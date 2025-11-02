function captureFrames(S, frames, filename)

    [~, name, ~] = fileparts(filename);
    videoName = [name, '_animation.mp4'];

    v = VideoWriter(videoName, 'MPEG-4');
    v.FrameRate = 30; % Set frame rate
    open(v);

    animateVelocityVectors(S, frames, v); 
end