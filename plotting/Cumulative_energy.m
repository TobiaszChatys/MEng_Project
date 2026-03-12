%% setup

clc; clear; close all;
data = load('Results/POD_data/POD_Results_L8_G6.mat');
energy = data.cumulative_energy;
%% Plotting
%--TODO: Add x and y labels as well as a legend

plot(1:length(energy), energy(:)' * 100, 'b-', 'Linewidth', 2)
hold on;

threshols =[0.90, 0.95, 0.99];
for threshold = 1:length(thresholds)
  
  yline(thresholds(threshold) * 100, '--');
  xline(modes_to_retain(threshold), ':');
  plot(modes_to_retain(threshold), thresholds(threshold) * 100, 'o');
  
end

