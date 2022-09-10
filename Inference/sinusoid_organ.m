clear;
%add audio path
addpath('Week3_audio');
[y, Fs] = audioread('organ.wav');
y = y(16000:20000)';

% Nength of data
N = length(y);
% Error variance
error_variance = 0.001;

%only retaining one side of information
signal_fft = fft(y);
P2 = abs(signal_fft/N);
P1 = P2(1:N/2+1); 
P1(2:end-1) = 2*P1(2:end-1);
P1 = log(P1);
f = Fs*(0:(N/2))/N;

figure(2)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of Organ.wav')
% xlim([0 22500])
xlabel('f (Hz)')
ylabel('Magnitude')


%find peaks and sort
[pks,locs] = findpeaks(P1, f,'MinPeakProminence', 3, 'MinPeakDistance',0.1);
%find frequency components of the peaks
fvalues = locs;
fvalues = 2* pi * fvalues / Fs;  


% Pass the identified frequencies into the model 
% Nength of data
N = length(y);
% Order of model
P = length(fvalues);
G = zeros(N, 2*P);
% Set bins of frequencies
w_bin = fvalues;


% Populate the model
for i =1:N
    for j = 1:P
        G(i, 2*j-1) = cos(i*w_bin(j));
        G(i, 2*j) = sin(i*w_bin(j));
    end
end

% ML estimation
theta_ML = inv(transpose(G)*G)*transpose(G)*y;
% MAP estimation
% Set theta prior parameters
% This is m_theta
theta_mean = zeros(2*P, 1);
% This is C_theta
theta_variance = eye(2*P);
[big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
    theta_mean,theta_variance,y);

% Estimate the output
y_ML = G*theta_ML;
y_MAP = G*theta_MAP;

% Plot results to see
figure(3)
plot(y);
hold on
plot(y_ML);
hold on
plot(y_MAP);
hold on 
xlabel('Time')
ylabel('Data')
legend('Original', 'ML', 'MAP')
title('Waveform of Original and Estimated Processes')
hold off

sound(y_MAP, Fs);
rootdirectory = "C:\Users\Youjing Yu\Desktop\SF1\Inference\Processed_audio";
file_name = fullfile(rootdirectory,'organ_MAP.wav');

audiowrite(file_name,y_MAP', Fs);






