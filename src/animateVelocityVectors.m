function animateVelocityVectors(S, frames, videoWriter, fig)

    if nargin < 4 || isempty(fig)
        fig = gcf;
    end
    
    isCapturing = exist('videoWriter', 'var') && ~isempty(videoWriter);
    
    % Initialize handles
    ax1 = []; ax2 = []; h = struct();

    for frame = 1:frames
        if ~isvalid(fig)
            break;
        end
        
        [X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);
        [Z1_masked, Z2_masked] = createMasks(X1, Y1, X2, Y2, Z1, Z2, X3, Y3);
        
        % Update graphics objects smoothly
        [ax1, ax2, h] = plotVectorMap(X1, Y1, U1, V1, X2, Y2, U2, V2, X3, Y3, Z1_masked, Z2_masked, ax1, ax2, h);
        
        title(ax2, sprintf('Velocity Vector Map - Frame %d', frame));
        
        if isCapturing
            frameCapture = getframe(fig);
            writeVideo(videoWriter, frameCapture);
        else
            drawnow limitrate;
            pause(0.01); 
        end
    end
end