function plotVectorMap(X1, Y1, U1, V1, X2, Y2, U2, V2, X3, Y3, Z1_masked, Z2_masked)

    % Create first axis air
    ax1 = axes("Position", [0.1, 0.1, 0.8, 0.8]);
    hold(ax1,'on')
    pcolor(ax1, X1 + 11, Y1, Z1_masked); 
    shading(ax1,'flat'); 
    colormap(ax1, hot)
    c1 = colorbar(ax1, "Position", [0.92, 0.1, 0.02, 0.8]);
    c1.Label.String = 'Air Velocity Magnitude (m/s)';

    % Create second axis liquid
    ax2 = axes('Position', [0.1, 0.1, 0.8, 0.8], 'Color', 'none');
    hold(ax2,'on')
    pcolor(ax2, X2 + 11, Y2, Z2_masked); 
    shading(ax2,'flat'); 
    colormap(ax2, jet);
    c2 = colorbar(ax2, "Position", [0.96, 0.1, 0.02, 0.8]);
    c2.Label.String = 'Liquid Velocity Magnitude (m/s)';

    linkaxes([ax1 ax2], 'xy');
    uistack(ax1, 'bottom')

 
    quiver(X1 + 11, Y1, U1, V1, 0.8, 'k');
    hold on;
    quiver(X2 + 11, Y2, U2, V2, 0.8, 'w');
    hold on; 
    plot(X3 + 11, Y3, 'g', 'LineWidth', 4)
    hold off;


    xlim([-15 15]);
    xlabel('X Position (mm)');
    ylabel('Y Position (mm)');
    title('L8, G3 Velocity Vector Map');
    ylim([0 28]);
    y_ticks = 0:2:28;
    yticks(y_ticks);
    yticklabels(y_ticks)
    xline(0, 'k--', 'LineWidth', 2);
end