clear; close all; clc; 

% Set the general settings:
audio = 'Audio/matin/Stim288.mat';
optimise = 0;
iterations = 1;
FS = 16e3;

% New implementation of the vocoder
[stim_CI, config, labels] = vocoder(audio, FS, 22, '22-electrodes',160, 'NOISE', 0, optimise, iterations);

% Save output
save('Audio/out/stim_CI_22Channels_160','FS','stim_CI','labels')