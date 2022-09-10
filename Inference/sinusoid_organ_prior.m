clear all;
%add audio path
addpath('Week3_audio');
[y, Fs] = audioread('organ.wav');
y = y(16000:20000)';

% Nength of data
N = length(y);
error_variance = 0.1;


%only retaining one side of information
signal_fft = fft(y);
P2 = abs(signal_fft/N);
P1 = P2(1:N/2+1); 
P1(2:end-1) = 2*P1(2:end-1);
P1 = log(P1);
f = Fs*(0:(N/2))/N;
figure(1)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of Organ.wav')
xlabel('f (Hz)')
ylabel('Magnitude')


%find peaks and sort
[pks,locs] = findpeaks(P1, f,'MinPeakProminence', 1.5, 'MinPeakDistance',0.1);
%find frequency components of the peaks
fvalues = locs;
fvalues = 2* pi * fvalues / Fs;  


% Pass the identified frequencies into the model 
% Length of  frequency data
N_freq = length(P1);
% Order of model; set to 1 since 1 frequency for each bin
P = 1;
% Set bins of frequencies
w_bin = fvalues;


% Incorporate prior so higher frequencies have lower amplitude
% Generate frequency bins; 1 bin for each frequency
num_bins = length(w_bin);
y_MAP_fre_list = [];
for i=1:num_bins
    signal = pks(i);
    theta_mean = 1/ signal;
    theta_variance = eye(1);
    G = 1;
    [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
    theta_mean,theta_variance, signal);
    % Estimate the output
    y_MAP = G*theta_MAP;    
    y_MAP_fre_list = [y_MAP_fre_list, y_MAP];
end


% Now we regenerate the original sound
G = zeros(N, 2*P);
% Populate the model
for i =1:N
    for j = 1:num_bins
        % Set prior on G according to the frequency values
        G(i, 2*j-1) = cos(i*w_bin(j))* y_MAP_fre_list(j);
        G(i, 2*j) = sin(i*w_bin(j))*y_MAP_fre_list(j);
    end
end


% MAP estimation
P = length(fvalues);
% Set theta prior parameters
% This is m_theta
theta_mean = zeros(2*P, 1);
% This is C_theta
theta_variance = eye(2*P);
[big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
    theta_mean,theta_variance,y);

% Estimate the output
y_MAP = G*theta_MAP;

% Plot results to see
figure(3)
plot(y);
hold on
plot(y_MAP);
hold on 
xlabel('Time')
ylabel('Data')
legend('Original', 'MAP')
title('Waveform of Original and Estimated Processes')
hold off

% sound(y, Fs);
sound(y_MAP,Fs);
rootdirectory = "C:\Users\Youjing Yu\Desktop\SF1\Inference\Processed_audio";
file_name = fullfile(rootdirectory,'organ_advanced_MAP.wav');

audiowrite(file_name,y_MAP', Fs);








