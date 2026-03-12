%% setup

clc; clear; close all;
data = load('Results/POD_data/POD_Results_L8_G6.mat');
energy = data.cumulative_energy;

%% color

color_liquid = [49, 116, 143] / 255; % #31748f

%% Plotting

figure,

plot(1:length(energy), energy(:)' * 100,'Color', color_liquid, 'Linewidth', 2)
hold on;

thresholds =[0.90, 0.95, 0.99];

for threshold = 1:length(thresholds)
  
  yline(thresholds(threshold) * 100, '--');
  xline(energy(threshold), '--');
  plot(energy(threshold), thresholds(threshold) * 100, 'o');
  
end

