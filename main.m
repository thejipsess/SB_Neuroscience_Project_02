clear; close all; clc; 

%Start of CI Model
% Set the general settings:
audio = 'Audio/in/Stim288.mat';
optimise = false;
iterations = 1;
FS = 16e3;
channels_list = [4, 8, 10, 16, 22];

% New implementation of the vocoder
[stim_CI, config, labels] = vocoder(audio, FS, 10, '10-electrodes',160, 'NOISE', 0, optimise, iterations);

% New implementation of the vocoder
for i = 1:length(channels_list)
    n_channels = channels_list(i);
    if n_channels > 9
        channel_spacing = strcat(num2str(n_channels),"-electrodes");
    else
        channel_spacing = 'default';
    end
    
    [stim_CI, config, labels] = vocoder(audio, FS, n_channels, channel_spacing,420, 'NOISE', 0, optimise, iterations);
    
    % Save output
    savename = strcat('Audio/out/numeric_stim_CI_', num2str(n_channels),'Channels_420');
    save(savename,'FS','stim_CI','labels');
    disp(strcat("saved: ", savename))
end

% Save output
save('Audio/out/stim_CI_10Channels_160.mat','FS','stim_CI','labels')


%Here data produced by CI model is loaded
features_file = 'Features/test_classifier.mat';
filename_in = 'Audio/out/stim_CI_10Channels_160.mat';
channels_list = [4, 8, 10, 16, 22];
start = 0;
finish = 96;

%this function uses matlab code, so index starts at 1
[~, features_spectral, labels] = load_data_optimized(filename_in, features_file, start+1, finish, FS);



%Here python code is runned
%Check if python is installed and change the location if needed
%pyenv('Version' , '/usr/bin/python3')
%When using Windows you can probably comment the next 2 lines
py.sys.setdlopenflags(int32(10));
py.importlib.import_module('ssl');

%Here the classifier is run. This uses python code
[scores, ~, ~] = classifier(features_spectral, labels, start, finish, 'rbf', labels(start+1:finish));
accuracy = mean(scores)*100;
boxplot(scores);