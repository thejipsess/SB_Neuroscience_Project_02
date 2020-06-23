clear
simOut = sim('CI_model2');

sim_time = input.Time;
in_sig1 = input.Data(1,:);
in_sig2 = input.Data(2,:);
in_sig3 = input.Data(3,:);
in_sig4 = input.Data(4,:);
in_sig5 = input.Data(5,:);
in_sig6 = input.Data(6,:);
in_sig7 = input.Data(7,:);
in_sig8 = input.Data(8,:);

Fs = size(sim_time, 1)/(sim_time(end));

full = in_sig1 + in_sig1 ++ in_sig1 + in_sig1 + in_sig1 + in_sig1 + in_sig1 + in_sig1;

%% Creat input

[y,Fs] = audioread('F:\Video edditing\Movie Voice KING\Movie Voice KING\Movie Voice mp3\Andrew\Ancient Mysteries\Do You Believe_ 1.mp3')
y = y(:,1);

audiowrite('sound',y, Fs);

save sound.mat y