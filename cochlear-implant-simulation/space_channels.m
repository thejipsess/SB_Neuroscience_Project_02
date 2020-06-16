function [Wn] = space_channels(n_channels, Fs, channel_spacing)

switch channel_spacing
    case 'linear'
        % Generate linearly spaced center frequencies
        center_freqs = transpose(repmat(linspace(357, 4740, n_channels), 2, 1));
        center_freqs = center_freqs/(Fs/2);
        
    case 'log'
        % Generate logarithmically spaced center frequencies
        center_freqs = transpose(repmat(logspace(log10(357), log10(4740), n_channels), 2, 1));
        center_freqs = center_freqs/(Fs/2);
        
        bandwidth = logspace(log10(357), log10(4740), n_channels)
        
    case '22-electrodes'
        % Create the cut-off frequencies for the bandpass filters
        Wn = [
        188, 313;
        313, 438;
        438, 563;
        563, 688;
        688, 813;
        813, 938;
        938, 1063;
        1063, 1188;
        1188, 1313;
        1313, 1568;
        1568, 1813;
        1813, 2063;
        2063, 2313;
        2313, 2588;
        2588, 3063;
        3063, 3563;
        3563, 4063;
        4063, 4688;
        4688, 5313;
        5313, 6063;
        6063, 6939;
        6939, 7938;
        ];
        
        % Normalise the cutoff frequencies  between 0 and 1, where 1 
        % corresponds to the Nyquist rate—half the sample rate
        Wn = Wn/(Fs/2);
        
    case 'default'
        if     n_channels==1
            Wn = repmat([4000], 1, 2);
            Bw=0.5*[-7500 7500]./(rate/2);
        elseif n_channels==2
            Wn = repmat([792; 3392], 1, 2);
            Bw=0.5*[-0984 0984; -4215 4215]./(rate/2);
        elseif n_channels==3
            Wn=repmat([0545; 1438; 3793], 1, 2);
            Bw=0.5*[-0491 0491; -1295 1295; -3414 3414]./(rate/2);
        elseif n_channels==4
            Wn=repmat([0460; 0953; 1971; 4078], 1, 2);
            Bw=0.5*[-0321 0321; -0664 0664; -1373 1373; -2842 2842]./(rate/2);
        elseif n_channels==5
            Wn=repmat([0418; 0748; 1339; 2396; 4287], 1, 2);
            Bw=0.5*[-0237 0237; -0423 0423; -0758 0758; -1356 1356; -2426 2426]./(rate/2);
        elseif n_channels==6
            Wn=repmat([0393; 0639; 1037; 1685; 2736; 4443],1, 2);
            Bw=0.5*[-0187 0187; -0304 0304; -0493 0493; -0801 0801; -1301 1301; -2113 2113]./(rate/2);
        elseif n_channels==7
            Wn=repmat([0377; 0572; 0866; 1312; 1988; 3013; 4565], 1, 2);
            Bw=0.5*[-0154 0154; -0234 0234; -0355 0355; -0538 0538; -0814 0814; -1234 1234; -1870 1870]./(rate/2);
        elseif n_channels==8
            Wn=repmat([0366; 0526; 0757; 1089; 1566; 2252; 3241; 4662], 1, 2);
            Bw=0.5*[-0131 0131; -0189 0189; -0272 0272; -0391 0391; -0563 0563; -0810 0810; -1165 1165; -1676 1676]./(rate/2);
        elseif n_channels==9
            Wn=repmat([0357; 0493; 0682; 0942; 1301; 1798; 2484; 3431; 4740], 1, 2);
            Bw=0.5*[-0114 0114; -0158 0158; -0218 0218; -0302 0302; -0417 0417; -0576 0576; -0796 0796; -1099 1099; -1519 1519]./(rate/2);
        else
            error('Wrong channel number for default channel_spacing');
        end
        
Wn = Wn/(rate/2);
Wn=Wn+Bw;
Wn(Wn>1) = 0.99;
Wn(Wn<0) = 0.01;
end
end