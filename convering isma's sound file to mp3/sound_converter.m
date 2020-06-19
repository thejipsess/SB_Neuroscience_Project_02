clear; close all; clc; 

% load in the work space
load('Stim288.mat')

% settings for the for loop
numfiles = length (labels);


for i=1:numfiles
    soundnumber= i
    sound = stim(soundnumber,:),fs
    
    name= labels{i};
    chr = mat2str(i)
    
    audiowrite(['audio/',name,chr,'.wav'], sound, fs)
end