clear; close all; clc; 

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
        [x,fs] = audioread(['Audio/in/',name,ext]);
    
        % Resample the audio
        x = resample(x, FS, fs);

        % Convert audio to the cochlear implant signal
        y = vocoder(x, FS, 22, 50, 'NOISE', 1);


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
    % load in the work space
    load('Audio/matin/Stim288.mat')

    % settings for the for loop
    numfiles = length (labels);
    stim_CI = zeros(length(labels),fs);
    
    
    for i=1:numfiles
        %selecting sounds
        soundnumber= i;
        sound = stim(soundnumber,:),fs;
    
        % Convert audio to the cochlear implant signal
        y = vocoder(sound, fs, 22, 50, 'NOISE', 0);
    
        % Plot the new audio signal commented out
        % figure;
        % plot((1:length(y))/fs, y);
        % xlabel('Time (second)');
        % title('Vocoded signal');

        % Normalise the sound signal between -1 and 1 since the audiowriter cannot
        % handle values outside this range.
        y2 = normalize(y, 'range', [-1 1]);
        
        % writing out the audio
        
        % name= labels{i};
        % chr = mat2str(i);
        % audiowrite(['Audio/matout/',name,chr,'.wav'], y2, fs)
        
        stim_CI(i,:) = y2;
        
        %counter
        count = i
    end
    
    save('Stim_CI','fs','stim_CI','labels','stim')
    
else 
    disp("error wrong input")    
end