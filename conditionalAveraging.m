[S, filename] = loadData('L8_G3.mat'); 
frame = 339;

%Get all data at once
[X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame);

% stats
h_max = max(Y3);
h_min = min(Y3);
h_mean = mean(Y3);
h_std = std(Y3);

fprintf('Film Height Statistics at Frame %d:\n', frame);
fprintf('Max Height: %.2f mm\n', h_max);
fprintf('Min Height: %.2f mm\n', h_min);
fprintf('Mean Height: %.2f mm\n', h_mean);
fprintf('Std Dev of Height: %.2f mm\n', h_std);

% Plot Velocity Vectors
quiver(X1 + 11, Y1, U1, V1, 0.8, 'k');
hold on;
quiver(X2 + 11, Y2, U2, V2, 0.8, 'b');
hold on;
plot(X3 + 11, Y3, 'g', 'LineWidth', 4)
hold off;
xlim([-15 15]);
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
title('L8-G3 Velocity Vector Map');
ylim([0 28]);
y_ticks = 0:2:28;
yticks(y_ticks);
yticklabels(y_ticks)
xline(0, 'k--', 'LineWidth', 2);