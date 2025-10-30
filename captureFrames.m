function captureFrames(S, frames)
    v = VideoWriter('velocity_vectors_animation.mp4', 'MPEG-4');
    v.FrameRate = 60; % Set frame rate
    open(v);

    animateVelocityVectors(S, frames, v); % Call the animation function
end