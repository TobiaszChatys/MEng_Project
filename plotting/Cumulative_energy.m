%% setup

clc; clear; close all;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
PROJ_ROOT  = fileparts(SCRIPT_DIR);
data = load(fullfile(PROJ_ROOT, 'scripts', 'Results', 'POD_data', 'POD_Results_L8_G6.mat'));
energy = data.cumulative_energy;

%% color

color_liquid = [49, 116, 143] / 255; % #31748f

%% Plotting

figure,

plot(1:length(energy), energy(:)' * 100,'Color', color_liquid, 'Linewidth', 2)
hold on;

thresholds = [0.90, 0.95, 0.99];

for threshold = 1:length(thresholds)

  modes_to_retain = find(energy >= thresholds(threshold), 1);
  yline(thresholds(threshold) * 100, '--');
  xline(modes_to_retain, ':');
  plot(modes_to_retain, thresholds(threshold) * 100, 'o');

end

xlabel('Mode Number');
ylabel('Cumulative Energy (%)');
