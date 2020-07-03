% set parameters
channels_list = [22, 16, 10, 8, 4];
Fs = 16e3;

spacing = nan(23,5);

for i = 1:length(channels_list)
    n_channels = channels_list(i);
    if n_channels > 9
        channel_spacing = strcat(num2str(n_channels),"-electrodes");
    else
        channel_spacing = 'default';
    end
    
    Wn = space_channels(n_channels, Fs, channel_spacing) * (Fs/2);
    
    boundaries = [Wn(1:end, 1); Wn(end,2)];
    spacing(1:n_channels+1, i) = [Wn(1,1); boundaries(2:end) - boundaries(1:end-1)];
end

channels_list_cat = categorical(channels_list);

for j = 1:length(channels_list)
    bp = barh(channels_list_cat(j), spacing(:,j), 'stacked');
    hold on
    bp(1).FaceColor = 'none';
    bp(1).EdgeColor = 'none';

    cmap = colormap(summer);
    for i = 2:length(spacing)
        coef = uint8(length(cmap)/channels_list(j));

        bp(i).FaceColor = cmap(coef*(i-1),:);
    end
end

ylabel('Number of channels');
xlabel('Frequency bands');
set(gca, 'FontSize', 18);


