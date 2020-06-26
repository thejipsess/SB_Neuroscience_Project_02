clear; close all; clc; 

%Start of CI Model
% Set the general settings:
audio = 'Audio/in/Stim288.mat';
optimise = false;
iterations = 1;
FS = 16e3;

n_channels = [4,8,10,16,22];
labels = repelem([0:9],200);
channel_spacing = {'default', 'default', '10-electrodes', '16-electrodes', '22-electrodes'};
files = {'numeric_stim_CI_4Channels_420', 'numeric_stim_CI_8Channels_420', 'numeric_stim_CI_10Channels_420', 'numeric_stim_CI_16Channels_420', 'numeric_stim_CI_22Channels_420'};
all_scores = zeros(5,250);
features = zeros(length(n_channels), 2000, 98, 4);
for i = 1:length(n_channels)
    % New implementation of the vocoder
    %[stim_CI, config, labels] = vocoder(audio, FS, n_channels(i), channel_spacing{i},160, 'NOISE', 0, optimise, iterations);

    % Save output
    %save('Audio/out/allOutputs_overwritten.mat','FS','stim_CI','labels')
    
    %Here data produced by CI model is loaded
    features_file = 'Features/allOutputs_overwritten.mat';
    filename_in = files{i};
    start = 0;
    finish = 2000;

    %this function uses matlab code, so index starts at 1
    [~, features_spectral, ~] = load_data_optimized(filename_in, features_file, start+1, finish, FS);
    features(i, :, :, :) = features_spectral;

    %Here python code is runned
    %Check if python is installed and change the location if needed
    %pyenv('Version' , '/usr/bin/python3')
    %When using Windows you can probably comment the next 2 lines
    py.sys.setdlopenflags(int32(10));
    py.importlib.import_module('ssl');

    %Here the classifier is run. This uses python code
    [scores, ~, ~] = classifier(features_spectral, labels, start, finish, 'rbf', labels(start+1:finish));
    all_scores(i,:) = scores;
end
%accuracy = mean(scores)*100;
boxplot(all_scores', n_channels);

%% 

