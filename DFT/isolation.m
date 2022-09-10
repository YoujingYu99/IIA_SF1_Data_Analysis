

%add audio path
addpath('Audio_Week1');

%vowel 'AE' in 'Captain' and consanant 'J' in 'John'
% [y1, Fs] = audioread('f1lcapae.wav');
%plot to see what the signal looks like
% plot(y1)
% xlabel('Time')
% ylabel('Amplitude')
% title('Time Response of Consonant')

% %now isolate the 'AE2
% y1_ae = y1(17000:17500, 1);
% plot(y1_ae);
% xlabel('Time')
% ylabel('Amplitude')
% title('Waveform of Vowel')
% %listen to make sure it is the vowel
% sound(y1_ae, Fs);
% now isolate the 'J'
% y1_j = y1(4000:5600, 1);
% plot(y1_j);
% xlabel('Time')
% ylabel('Amplitude')
% title('Waveform of Consonant')
%listen to make sure it is the consonant
% sound(y1_j, Fs);
%plot(y1_j);


%steady and transient notes in grosse
[y2, Fs] = audioread('piano_clean.wav')
%plot to see what the signal looks like
% plot(y2)
% xlabel('Time')
% ylabel('Amplitude')
% title('Time Response of Consonant')

% %now isolate the steady note
y2_steady = y2(12100:12300,1);
% plot(y2_steady)
% xlabel('Time')
% ylabel('Amplitude')
% title('Waveform of Steady Note')
% %listen to make sure it is the steady note
% sound(y2_steady, Fs);
% %now isolate the transient note
y2_tran = y2(8800:9300,1);
plot(y2_tran)
xlabel('Time')
ylabel('Amplitude')
title('Waveform of Transient Note')
% %listen to make sure it is the transient note
% sound(y2_tran, Fs);



% %try out different windows
% %choose hamming window with length N = N
% Nsamps = length(y1_ae);
% audio_window(y1_ae, 'hamm', Nsamps, Fs);











