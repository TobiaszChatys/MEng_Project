function animateVelocityVectors(S, frames, videoWriter, fig)

    isCapturing = exist('videoWriter', 'var');

    for frame = 1:frames

        fig = figure(fig);
        [X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);
        
        quiver(X1 + 11, Y1, U1, V1, 0.8, 'k');
        hold on;
        quiver(X2 + 11, Y2, U2, V2, 0.8, 'b');
        hold on; 
        plot(X3 + 11, Y3, 'g', 'LineWidth', 4)
        hold off;

        xlabel('X Position (mm)');
        ylabel('Y Position (mm)');
        title(sprintf('L8-G3 Velocity Vector Map - Frame %d', frame));
        
        xlim([-15 15]);
        ylim([0 28]);

        y_ticks = 0:2:28;
        yticks(y_ticks);
        yticklabels(y_ticks)
        xline(0, 'k--', 'LineWidth', 2);

        if isCapturing
            frameCapture = getframe(fig);
            writeVideo(videoWriter, frameCapture);
        else
            pause(1/60);  % 60 fps
        end
    end
end