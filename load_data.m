function [] = load_data(filename)
load('Stim288.mat');

range = 1:96;

EE1_list = zeros(length(range),98,16001);
EE2_list = zeros(length(range),98,16001);
EE3_list = zeros(length(range),98,16001);
EE4_list = zeros(length(range),98,16001);

%for i = 1:size(stim,1)
for i = range
    i
    s = stim(i,:)';
    AC_Model;
    EE1_list(i,:,:) = EE1; 
    EE2_list(i,:,:) = EE2;
    EE3_list(i,:,:) = EE3; 
    EE4_list(i,:,:) = EE4; 
end
 save(filename, 'EE1_list', 'EE2_list', 'EE3_list', 'EE4_list');
end

