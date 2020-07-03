function [CI_signals, cf_rand, labels] = vocoder(x, Fs, n_channels, channel_spacing, cutoff, vocoder_type, verbose, optimise, iterations)
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

if ischar(x)
    % Check if the audio input is specified as a .mat files
    if strcmp(x(end-3:end), '.mat')
        % Load the .mat file
        load(x)
        
        % Transpose the audio matrix
        audio = transpose(stim);
        
        % Initialise cell array x
        x = cell(1);
        
        % Convert audio to a cell array
        for i = 1:size(audio, 2)
            x{i} = audio(:,i);
        end
        
        disp('loaded .mat file succesfully')
        % Still need to incorporate and else statement to handle individual
        % files or folders filled with audio files.
    elseif isfolder(x)
        % Load all audio files from folder
        folder = x;
        if folder(end) ~= '/'
            folder = [folder, '/'];
        end
        
        % Set supported files types
        filetypes = {'.mp3', '.wav', '.flac'};
        
        % Count number files in folder
        files = dir(folder);
        numfiles = length(files);
        
        % Count number of supported files in folder
        numfiles_sup = length(dir(strcat(folder, '*.mp3')));
        numfiles_sup = numfiles_sup + length(dir(strcat(folder, '*.wav')));
        numfiles_sup = numfiles_sup + length(dir(strcat(folder, '*.flac')));
        
        % Initliase x to store all audio files
        x = cell(1);
        
        index = 1;
        for i=1:numfiles

            % Set settings
            FS = 16e3;
            % filetype = '.mp3';
            % filename = 'beethoven_fur_elise_orig';

            % getting file name withouth extentsion
            file = files(i,1);
            file = file.name;
            [filepath,name,ext] = fileparts(file);
            
            % Check if the filetype is supported
            if ismember(ext, filetypes)
                % Load  the audio file
                [audio,fs] = audioread([folder,name,ext]);

                % Resample the audio
                audio = resample(audio, FS, fs);
                
                % Ensure the audio is mono
                audio = mean(audio, 2);
                
                % Store the audio in 
                x{index} = audio;
                index = index + 1;
            end
        end
        disp(strcat("Loaded ", num2str(numfiles_sup), ' files.'))
        labels = 'no labels available';
    else
        labels = 'no audio could be loaded';
    end
end

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
else
    cf_rand = 'not optimised';
    iterations = 1;
end

% Initialise the CI_signals variable which will store all the different
% channel configurations that are being used, and with them stores all the
% CI converted audio signals.
CI_signals = cell(iterations, 2);
for arch = 1:iterations
    CI_signals{arch,2} = cf_rand(:,arch);
end

% Loop though all audio segments in x (as = Audio Segment)
for as = 1:length(x)
    % Check length of audio signal
    npts = length(x{as});
    
    % Initialise y, which will contain the processed CI audio signal for
    % all channel stetups
    y = NaN(length(x{as}), iterations);

    % Apply high-pass pre-emphasis filter
    pre=0.9378;
    xx=filter([1 -pre], 1, x{as})';

    index = 1;
    for iter = 1:iterations
        if optimise
            % Select current channel setup
            Wn = Wn_rand(:,index:index+1);
        end
        
%         % Show butterworth filter design
%         if verbose
%             freqz(blp, alp)
%         end
        
        % Generate noise carrier (only for noise vocoders)
        if( strcmp(vocoder_type, 'NOISE') )
            noise = rand( length(x{as}),1 );
            noise = noise(:);
        end

        % Initialse vocoded_x
        vocoded_x=zeros(npts,1);

        for i=1:n_channels

            % Find the filter coefficients for each bandpass filter
            [b,a] = butter(4,Wn(i,:));
            
            % Show butterworth filter design
            if verbose
                figure
                freqz(b, a)
            end

            % now filter the input waveform using filter #i
            filtwav = filtfilt(b,a,xx)';

            % Half-wave rectification
            filtwav(filtwav<0) = 0;
            
            % GeneFs lowpass filter coefficients (for envelope extraction):
            fc = cutoff /(Fs/2);
            % The butter-worth filter create a frequency response as flat as possible
            % in the pass band region. Here it is a low pass filter.
            [blp,alp]=butter(2,fc,'low'); % geneFs filter coefficients

            % Filter the band-pass filtered signal to extract its envelope
            % (Overall shape)
            envelope = filter(blp,alp,filtwav);
            envelope = envelope(:);

            % If noise vocoder is selected, then multiply the envelope with the
            % noise carrier.
            % Basically, we modulate the noise by the envelope
            if(strcmp(vocoder_type, 'NOISE'))
                % Ensure that the envelope and noise are of the same length
                if length(envelope) > length(noise)
                    envelope(length(noise)+1:length(envelope))=[];
                elseif length(noise) > length(envelope)
                    noise(length(envelope)+1:length(noise))=[];
                end
                
                % Normalise noise signal
                source = noise./(max(abs(noise)));
                
                % Amplitude modulate the noise carrier with the envelope
                % and perform the bandpass filter
                fn = filtfilt(b,a,envelope.*source);

            else
                %    If tone vocoder is selected, then multiply the envelope with a
                %    tone carrier.
                %    Basically, we modulate the tone by the envelope

                %     Tone with freq. at the center of the band-pass filter
                f = exp(mean(log(Wn(i,:)))) * (Fs/2);
                tone = sin(2*pi*(1:length(envelope))*f/Fs)';
                tone = tone(:);
                fn = envelope.*tone;
            end

            % sum bands with equal gain in each channel
            vocoded_x = vocoded_x + fn;
            
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

        end
        y(:, iter) = vocoded_x;
        index = index + 2;
    end
        
    for i = 1:size(y, 2)
        % Scale output waveform to have same rms as original
        y(:,i) = y(:,i) * (rms(x{as})/rms(y(:,i)));

        % Normalise the sound signal such that it does not surpass the
        % extreme valuse of the input sound. This ensure that the
        % loudness of both signals is comparable.
        y(:,i) = normalize(y(:,i), 'range', [min(x{as}) max(x{as})]);
        
        % Store the vocoded audio signal in the CI_signals variable
        CI_signals{i,1}{1, as} = y(:,1);
    end
end

% Stop timer
toc

end