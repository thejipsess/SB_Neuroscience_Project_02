%% Memory inefficient version, might delete later

filename = 'AC_output_voice.mat';
%
load_data(filename);
load(filename);

[temporal_mod_EE1, spectral_mod_EE1] = feature_extraction(EE1_list);
[temporal_mod_EE2, spectral_mod_EE2] = feature_extraction(EE2_list);
[temporal_mod_EE3, spectral_mod_EE3] = feature_extraction(EE3_list);
[temporal_mod_EE4, spectral_mod_EE4] = feature_extraction(EE4_list);

feature_file = 'features_voice_output.mat';
save(feature_file, 'freq_over_mfr_EE1', 'freq_over_mfr_EE2', 'freq_over_mfr_EE3', 'freq_over_mfr_EE4', ...
                   'mfr_over_time_EE1', 'mfr_over_time_EE2', 'mfr_over_time_EE3', 'mfr_over_time_EE4');

%% Optimized version to safe memory
%This version does not load all the data but loads the data per sound clip
%and only stores the extracted features.

features_file = 'test_classifier4.mat';
input_filename = 'stim_CI_4Channels.mat';
start = 0;
finish = 3; %96

%this function uses matlab code, so index starts at 1
[~, features_spectral, labels] = load_data_optimized(input_filename, features_file, start+1, finish);
load(features_file, 'features_spectral');
load('stim_CI_4Channels.mat')

%Here python code is runned
%Check if python is installed and change the location if needed
%pyenv('Version' , '/usr/bin/python3')
py.sys.setdlopenflags(int32(10));
py.importlib.import_module('ssl');

scores = classifier(features_spectral,labels, start, finish, 'rbf');

%% To load wav files

%[y, Fs] = audioread('0_jackson_0.wav')