function space_channels(n_channels, Fs, type)

if strcmp(type, 'linear')
center_freqs = transpose(repmat(linspace(357, 4740, n_channels), 2, 1));
center_freqs = center_freqs/(Fs/2);

elseif strcmp(type, 'log')
center_freqs = transpose(repmat(logspace(log10(357), log10(4740), n_channels), 2, 1));
center_freqs = center_freqs/(Fs/2);
end

Wn=repmat([0357; 0493; 0682; 0942; 1301; 1798; 2484; 3431; 4740], 1, 2);

end