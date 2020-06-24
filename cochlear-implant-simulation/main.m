clear; close all; clc; 

% Set the general settings:
audio = 'Audio/in/Stim288.mat';
optimise = 1;
iterations = 3;
FS = 16e3;

% New implementation of the vocoder
[stim_CI, config, labels] = vocoder(audio, FS, 10, '10-electrodes',160, 'NOISE', 0, optimise, iterations);

% Save output
save('Audio/out/stim_CI_10Channels_160','FS','stim_CI','labels')