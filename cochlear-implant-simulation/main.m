clear; close all; clc; 

% Set settings
FS = 16e3;
filetype = '.mp3';
filename = 'Japanese jo';

% Load  the audio file
[x,fs] = audioread(['Audio/in/',filename, filetype]);

% Resample the audio
x = resample(x, FS, fs);

% Convert audio to the cochlear implant signal
y = vocoder(x, FS, 22, 50, 'NOISE', 1);


% play sound
sound(y, FS)

% Plot the new audio signal
figure;
plot((1:length(y))/FS, y);
xlabel('Time (second)')
title('Vocoded signal')

% Normalise the sound signal between -1 and 1 since the audiowriter cannot
% handle values outside this range.
y2 = normalize(y, 'range', [-1 1]);

% Export the converted audio file
audiowrite(['Audio/out/',filename,'.wav'], y2, FS);