%% Video Export Script for Presentation
% This script exports the Velocity Vector Map animation as a 16:9 MPEG-4 video.
% It includes high-resolution heatmaps, quivers, and the "Crap Vector" callout.

clear; clc; close all;
addpath('src/');
addpath('data/Cases/');

% 1. Load the dataset
fprintf('Loading data...\n');
[S, filename] = loadData('L8_G9.mat');

% 2. Configuration
numFrames = 300; % Set the number of snapshots to export (e.g., 2,500 for full run)
fprintf('Starting video export for %d frames...\n', numFrames);

% 3. Capture and Export
% This calls captureFrames which sets 1920x1080 resolution and MPEG-4 format
captureFrames(S, numFrames, filename);

fprintf('Done! Check the project folder for the .mp4 file.\n');
