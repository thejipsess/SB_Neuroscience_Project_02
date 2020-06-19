function [vocoded_x, cf_rand] = vocoder(x, Fs, n_channels, channel_spacing, cutoff, vocoder_type, verbose, optimise, iterations)
% Implementation of cochlear implant simulation, referred to as vocoder.
% The first argument is the signal.
% The second argument is the sampling Fs (preferred 16Khz)
% The third argument is the number of spectral channels between 2 and 9.
% The fourth argument specifies the type of channel spacing
% The fifth argument is the cutoff frequency for envelope extraction. The higher the
%   cutoff, the more fine structure will be available in the vocoded signal.
% The sixth argument specifies the type of vocoder: either NOISE or TONE vocoders.
% The seventh argument provides more details regarding the bandpass
% filters.
% The 8th argument states whether channel placement should be optimised
% The 9th argument states how many iterations should be performed for
% optimisation
% See the example file 'example.m' on how to use this code.
% Refer to Shannon, Zeng, Kamath, Wygonski and Ekelid (1995). 
% Speech Recognition with Primarily Temporal Cues, Science
% for more information about vocoders. 
% 
% Code written by Vahid Montazeri, October 2018.
% Code editted by Jip de Kok, June 2020.

% Start timer
tic

% If not all inputs are defined then set missing inputs to their defaults.
switch nargin
    case 2
        n_channels = 8;
        channel_spacing = 'default';
        cutoff = 160;
        vocoder_type = 'NOISE';
        verbose = 1;
        optimise = 0;
        iterations = 1;
    case 3
        channel_spacing = 'default';
        cutoff = 160;
        vocoder_type = 'NOISE';
        verbose = 1;
        optimise = 0;
        iterations = 1;
    case 4
        cutoff = 160;
        vocoder_type = 'NOISE';
        verbose = 1;
        optimise = 0;
        iterations = 1;
    case 5
        vocoder_type = 'NOISE';
        verbose = 1;
        optimise = 0;
        iterations = 1;
    case 6
        verbose = 1;
        optimise = 0;
        iterations = 1;
    case 7
        optimise = 0;
        iterations = 1;
    case 8
        if optimise
            iterations = 50;
        else
            iterations = 1;
        end
end

% Make sure the signal is a column vector
if( size(x,2)>size(x,1) )
    x = x';
end
x = x(:,1);

npts=length(x);

% Apply high-pass pre-emphasis filter
pre=0.9378;
xx=filter([1 -pre], 1, x)';

% Center Freq. based on number of channels. These values are taken from
% Dorman, Loizou, Rainey, (1997). Speech intelligibility as a function of the number of channels
% of stimulation for signal processors using sine-wave
% and noise-band outputs. JASA

Wn = space_channels(n_channels, Fs, channel_spacing);

if optimise
    verbose = 0;

end

% Randomise channel placement for optimisation
if optimise
    [Wn_rand, cf_rand] = randomise_channels(Wn, iterations);
    vocoded_x = NaN(length(x), iterations);
else
    cf_rand = 'not optimised';
end

index = 1;
for iter = 1:iterations
    Wn = Wn_rand(:,index:index+1);
    % GeneFs lowpass filter coefficients (for envelope extraction):
     fc=cutoff /(Fs/2);
     % The butter-worth filter create a frequency response as flat as possible
     % in the pass band region. 
    [blp,alp]=butter(2,fc,'low'); % geneFs filter coefficients

    % GeneFs noise carrier only for noise vocoders
    if( strcmp(vocoder_type, 'NOISE') )
        noise = rand( length(x),1 );
        noise = noise(:);
    end

    vocoded_x=zeros(npts,1);

    for i=1:n_channels

        %     Find the filter coefficients for each bandpass filter
        [b,a] = butter(4,Wn(i,:));

        if(verbose)
            [h,f]=freqz(b,a,1024,Fs);        
            figure(1);
            plot(f,20*log10(abs(h+eps)),'LineWidth',2);
            hold on;
            min_f = min(Wn(:,1))*Fs/2 - 50;
            axis([min_f Fs/2 -6 0.5]);
            x_ticks = round(logspace(log10(min_f), log10(Fs/2), n_channels));
            set(gca,'XScale','log', 'XTick', x_ticks);                  
            xlabel('Frequency (Hz)');
            ylabel('Filter gain (dB)');
            title('Band-pass filters transfer function');
            grid on;
        end

        % now filter the input waveform using filter #i
        filtwav = filtfilt(b,a,xx)';

        %     Half-wave rectification
        filtwav(filtwav<0) = 0;

        %     Filter the band-passed filtered signal to extract its envelope
        %     (Overall shape)
        envelope=filter(blp,alp,filtwav);
        envelope = envelope(:);

        %     If noise vocoder is selected, then multiply the envelope with the
        %     noise carrier.
        %    Basically, we modulate the noise by the envelope
        if(strcmp(vocoder_type, 'NOISE'))
            if length(envelope) > length(noise)
                envelope(length(noise)+1:length(envelope))=[];
            end

            if length(noise) > length(envelope)
                noise(length(envelope)+1:length(noise))=[];
            end

            source = noise./(max(abs(noise)));
            fn=filtfilt(b,a,envelope.*source);

        else
            %     If tone vocoder is selected, then multiply the envelope with a
            %     tone carrier.
            %    Basically, we modulate the tone by the envelope

            %     Tone with freq. at the center of the band-pass filter
            f = exp(mean(log(Wn(i,:)))) * (Fs/2);
            tone=sin(2*pi*(1:length(envelope))*f/Fs)';
            tone = tone(:);
            fn = envelope.*tone;

        end

        % sum bands with equal gain in each channel
        vocoded_x = vocoded_x + fn;
        
    end
    vocoded_x(:, iter) = vocoded_x;
    index = index + 2;
end
    % Scale output waveform to have same rms as original
    vocoded_x = vocoded_x * (rms(x)/rms(vocoded_x));
    
    % Stop timer
    toc
end