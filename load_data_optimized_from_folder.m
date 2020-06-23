function [features_spectral, labels] = load_data_optimized_from_folder(filename,x)
%load('Stim288.mat');

files = dir(sprintf('recordings/%d_*.wav', x));

range = 1:size(files, 1);
channels = 98;
labels = repelem(x, size(files,1));

features_spectral = zeros(length(range), channels, 4);
for i = range
    fprintf('digit %d, file %d/%d\n',x,  i, size(files, 1));
    [s, fs] = audioread(join([files(i).folder, "/", files(i).name], ""));
    AC_Model;
    [~, features_spectral(i,:,1)] = feature_extraction(reshape(EE1, 1,channels,[]));
    [~, features_spectral(i,:,2)] = feature_extraction(reshape(EE2, 1,channels,[]));
    [~, features_spectral(i,:,3)] = feature_extraction(reshape(EE3, 1,channels,[]));
    [~, features_spectral(i,:,4)] = feature_extraction(reshape(EE4, 1,channels,[]));
end
 save(filename, 'features_spectral');
end

