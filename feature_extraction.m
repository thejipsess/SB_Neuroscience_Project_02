function [temporal_mod, spectral_mod] = feature_extraction(EE_list)

temporal_mod = zeros(size(EE_list,1),16001);
for i = 1:size(EE_list,1)
    temporal_mod(i,:) = mean(squeeze(EE_list(i, :, :)),1);
end

spectral_mod = zeros(size(EE_list,1),98);
for i = 1:size(EE_list,1)
    spectral_mod(i,:) = mean(squeeze(EE_list(i, :, :)),2);
end

end
