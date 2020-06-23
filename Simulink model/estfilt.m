function estfilt(nChannels,bpfcutoff)

% Function: estfilt.m
%
% estimate the parameters of the bandpass filters--
%
% Advanced Bionics, Clarion, Platinum Sound Processor Filter Table (2000) CIS
%
% input: nChannels - number of channels
%
% adapted from Loizou's estfilt.m
%
% Rebecca Reich
% January 2002


global filterA filterB center Srate
SAVE=0; % if 1, save center frequencies and bandwidths in a file
FS=Srate/2;
nOrd=6;
center=zeros(1,nChannels);
% ========== frequency boundaries for bandpass filters, depending on # of channels ========== %
if nChannels == 16
% ============ 16-channels, logarithmically-spaced ======================
UpperFreq=6800; LowFreq=bpfcutoff;
range=log10(UpperFreq/LowFreq);
interval=range/nChannels;
center=zeros(1,nChannels);
for i=1:nChannels % ----- Figure out the center frequencies for all channels
upper1(i)=LowFreq*10^(interval*i);
lower1(i)=LowFreq*10^(interval*(i-1));
end
else
switch nChannels
case 1
upper1 = 6800;
lower1 = bpfcutoff;
case 2
upper1 = [1387 6800];
lower1 = [bpfcutoff 1387];
case 3
upper1 = [877 2196 6800];
lower1 = [bpfcutoff 877 2196];
case 4
upper1 = [697 1387 2762 6800];
lower1 = [bpfcutoff 697 1387 2762];
case 5
upper1 = [607 1053 1827 3170 6800];
lower1 = [bpfcutoff 607 1053 1827 3170];
case 6
upper1 = [554 877 1387 2196 3475 6800];
lower1 = [bpfcutoff 554 877 1387 2196 3475];
case 7
upper1 = [519 769 1140 1689 2504 3711 6800];
lower1 = [bpfcutoff 519 769 1140 1689 2504 3711];
case 8
upper1 = [494 697 983 1387 1958 2762 3898 6800];
lower1 = [bpfcutoff 494 697 983 1387 1958 2762 3898];
end
end
for i = 1:nChannels % take GEOMETRIC MEAN (not simple average)
center(i)=sqrt(upper1(i)*lower1(i));
end
% ========== for SAVE = 1 print filter values to a file ========== %
if SAVE==1
fps=fopen('filters.txt','a+');
fprintf(fps,'%d channels:\n',nChannels);
for i=1:nChannels
fprintf(fps,'%d ',round(upper1(i)-lower1(i))); %bandwidths
end

fprintf(fps,'\n');
for i=1:nChannels
fprintf(fps,'%d ',round(center(i))); %center frequencies
end
fprintf(fps,'\n=======\n');
fclose(fps);
end
% ========== check for aliasing ========== %
if FS<upper1(nChannels) % need sRate >= 2 * max freq
useHigh=1;
else
useHigh=0;
end
% ========== design the filters ========== %
filterA=zeros(nChannels,nOrd+1);
filterB=zeros(nChannels,nOrd+1);
for i=1:nChannels
W1=[lower1(i)/FS, upper1(i)/FS];
if i==nChannels
if useHigh==0
[b,a]=butter(0.5*nOrd,W1);
else
[b,a]=butter(nOrd,W1(1),'high');
end
else
[b,a]=butter(0.5*nOrd,W1);
end
filterB(i,1:nOrd+1)=b;
filterA(i,1:nOrd+1)=a;
end

end