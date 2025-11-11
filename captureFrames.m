function captureFrames(S, frames, filename)

    [~, name, ~] = fileparts(filename);
    videoName = [name, '_animation.mp4'];

    % high resolution figure:
    fig = figure("Position",[100 100 1920 1080],"color","w","Visible","off");


    v = VideoWriter(videoName, 'MPEG-4');
    v.FrameRate = 30; % Set frame rate
    v.Quality = 100;
    open(v);

    animateVelocityVectors(S, frames, v, fig); 
    close(v);
    close(fig);

end