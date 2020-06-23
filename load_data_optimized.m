function [features_temporal, features_spectral, labels] = load_data_optimized(filename_in, filename_out, start, finish, FS)
%load('Stim288.mat');
load(filename_in);

range = start:finish;
channels = 98;
T = 16001;
stim = stim_CI{1,1}';
fs=FS;

features_temporal = zeros(length(range), T, 4); %4 since we want to store EE1 to EE4
features_spectral = zeros(length(range), channels, 4);
for i = range
    i
    s = stim(i,:)';
    AC_Model;
    [features_temporal(i,:,1), features_spectral(i,:,1)] = feature_extraction(reshape(EE1, [1,channels,T]));
    [features_temporal(i,:,2), features_spectral(i,:,2)] = feature_extraction(reshape(EE2, [1,channels,T]));
    [features_temporal(i,:,3), features_spectral(i,:,3)] = feature_extraction(reshape(EE3, [1,channels,T]));
    [features_temporal(i,:,4), features_spectral(i,:,4)] = feature_extraction(reshape(EE4, [1,channels,T]));
end
 save(filename_out, 'features_temporal', 'features_spectral');
end

