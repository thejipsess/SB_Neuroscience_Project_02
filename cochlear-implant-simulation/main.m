clear; close all; clc; 

% Set the general settings:
audio = 'Audio/matin/Stim288.mat';
optimise = 0;
iterations = 1;



prompt = 'Do you have a .mat file? y/n:';
Userinput = input(prompt,'s');
USerinput = convertCharsToStrings(Userinput)

if strcmp(Userinput,'n')

    % setting audio in map
    directory = dir('Audio/in/*.mp3');
    numfiles = length (directory);
    

    for i=1:numfiles

        % Set settings
        FS = 16e3;
        % filetype = '.mp3';
        % filename = 'beethoven_fur_elise_orig';
    
        % getting file name withouth extentsion
        file = directory(i,1);
        file = file.name;
        [filepath,name,ext] = fileparts(file);

        % Load  the audio file
        [sound,fs] = audioread(['Audio/in/',name,ext]);
    
        % Resample the audio
        x = resample(x, FS, fs);

        % Convert audio to the cochlear implant signal
        y = vocoder(x, FS, 22, '22-electrodes', 240, 'NOISE', 1);
        
        % New implementation of the vocoder (DOESN'T WORK YET WITH OPTIMISATION)
        %[y, config] = vocoder(sound, FS, 22, '22-electrodes', 240, 'NOISE', 0, 1, 100);


        % play sound commented out
        % sound(y, FS)

        % Plot the new audio signal
        figure;
        plot((1:length(y))/FS, y);
        xlabel('Time (second)');
        title('Vocoded signal');

        % Normalise the sound signal between -1 and 1 since the audiowriter cannot
        % handle values outside this range.
        y2 = normalize(y, 'range', [-1 1]);

        % Export the converted audio file
        filename = sprintf(name)
        audiowrite(['Audio/out/',filename,'.wav'], y2, FS);
        
        %counter
        count = i
    end
    
    
elseif strcmp(Userinput,'y')
    % Set settings
    FS = 16e3;

    % New implementation of the vocoder
    [stim_CI, config] = vocoder(audio, FS, 22, '22-electrodes', 240, 'NOISE', 0, optimise, iterations);

    % writing out the audio

    % name= labels{i};
    % chr = mat2str(i);
    % audiowrite(['Audio/matout/',name,chr,'.wav'], y2, fs)

    save('Audio/out/stim_CI_22Channels','fs','stim_CI','labels')
    
else 
    disp("error wrong input")    
end

