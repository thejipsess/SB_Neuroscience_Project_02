clear; close all; clc; 

% Set settings
FS = 16e3;
filetype = '.mp3';
filename = 'Together They Found the Courage 2';
filename = 'beethoven_fur_elise_orig';
filename = 'sail';

% Load  the audio file
[x,fs] = audioread(['Audio/in/',filename, filetype]);

% Resample the audio
x = resample(x, FS, fs);

% Convert audio to the cochlear implant signal
y = vocoder(x, FS, 4, 50, 'NOISE', 1);
y2 = vocoder(x, FS, 4, 4100, 'NOISE', 1);
y = vocoder(x, FS, 22, '22-electrodes', 160, 'NOISE', 1);

% Plot the new audio signal
figure;
plot((1:length(x))/FS, x(:,1));
hold on
plot((1:length(y))/FS, y);
xlabel('Time (second)')
title('Vocoded signal')

% Normalise the sound signal between -1 and 1 since the audiowriter cannot
% handle values outside this range.
y2 = normalize(y, 'range', [-1 1]);

% Export the converted audio file
audiowrite(['Audio/out/',filename,'.wav'], y2, FS);