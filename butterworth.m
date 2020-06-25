[d,s] = xlsread('data_excel1.xls');
t = datenum(s(2:end,1), 'yyyy-mm-dd HH:MM:SS');
figure(1)
plot(t, d)
grid
datetick('x', 'HH:MM:SS', 'keepticks')
Ts = mean(diff(t))*(24*60*60);                          % Sampling Time (sec)
Fs = 1/Ts;                                              % Sampling Frequency (Hz)
Fn = Fs/2;                                              % Nyquist Frequency (Hz)
L = length(t);
dsm = d-mean(d);                                        % Subtract Mean
FTd = fft(dsm)/L;                                       % Fourier Transform
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                     % Frequency Vector
Iv = 1:length(Fv);                                      % Index Vector
figure(2)
plot(Fv, 2*abs(FTd(Iv,:)))
grid
axis([0  0.01    ylim])
Wp = [0.0002  0.0070]/Fn;                               % Specify Bandpass Filter
Ws = [0.0001  0.0075]/Fn;                               % Stopband (normalised)
Rp = 10;                                                % Passband Ripple (dB)
Rs = 30;                                                % Stopband Ripple (dB)
[n,Wn] = buttord(Wp, Ws, Rp, Rs);
[z,p,k] = butter(n,Wn);                                 % ZPK 
[SOS,G] = zp2sos(z,p,k);                                % Convert To SOS for Stability
figure(3)
freqz(SOS, 2^16, Fs)                                    % Filter Bode Plot
d_filt = filtfilt(SOS, G, d);                           % Filter Signal
figure(4)
plot(t, d_filt)
grid
datetick('x', 'HH:MM:SS', 'keepticks')