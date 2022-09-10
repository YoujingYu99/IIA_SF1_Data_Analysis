clear;
%add audio path
addpath('Week3_audio');
[y_voice, Fs1] = audioread('f1lcapae.wav');

%now isolate the AE2
y_ae = y_voice(17000:17500, 1);
power_ae = log(abs(fft(y_ae, Fs1)).^2);
y_j = y_voice(4000:5600, 1);
power_j = log(abs(fft(y_j, Fs1)).^2);


[y_note, Fs2] = audioread('piano_clean.wav');
y_note_steady = y_note(12100:12300,1);
y_note_tran = y_note(8800:9300,1);

% Select the audio to analyse
y_original = y_j;
error_variance = 0.1;
N = length(y_original);


ML_MSE_list = [];
MAP_MSE_list = [];
model_LLR_list = [];

% MSE for different model orders
for i = 1:200
    % Set different AR orders
    G_AR = ARmodel(y_original, N,i);
    G = G_AR;
    y = y_original(i+1:end);
    % ML estimation
    theta_ML = inv(transpose(G)*G)*transpose(G)*y;
    y_ML = G*theta_ML;
    ML_MSE = mean((y_ML - y).^2);
    ML_MSE_list = [ML_MSE_list, ML_MSE];

    %MAP estimation
    % This is m_theta; column vector of length i(ith order)
    theta_mean = zeros(i,1);
    % This is C_theta; identity vector of size i
    theta_variance = 10*eye(i);
    % Posterior parameters
    % Use the previous error values
    if i >1
      error_variance = MAP_MSE_list(i-1);
    else
      error_variance = error_variance;
    end
    [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
        theta_mean,theta_variance,y);
    y_MAP = G*theta_MAP;
    MAP_MSE = mean((y_MAP - y).^2);
    MAP_MSE_list = [MAP_MSE_list, MAP_MSE];

    % Model Selection
    marginal_LLR = model_LLR(N,i,G,error_variance, theta_mean,theta_variance,y);
    model_LLR_list = [model_LLR_list; marginal_LLR];


end

% % Clean up LLR list by removing infinity and NaN values
% model_LLR_list(isinf(model_LLR_list)|isnan(model_LLR_list)) = 0;

% Return the optimum offset and the LLR
[optimum_LLR, optimum_order] = max(model_LLR_list);
% Convert back from the log domain to normal and normalise
posterior_probs = exp(model_LLR_list);
posterior_probs_norm = posterior_probs /sum(posterior_probs);

% Plot posterior probabilities against positions
figure
n = linspace(1, 100, 100);
plot(n, model_LLR_list)
xlabel('Model Order')
ylabel('Probability Density')
title('Posterior Probabilities Against Order')

% Determine optimum order
P = optimum_order;
G_AR = ARmodel(y_original, N,P);
G = G_AR;
y = y_original(P+1:end);
% ML estimation
theta_ML = inv(transpose(G)*G)*transpose(G)*y;

%MAP estimation
% This is m_theta; column vector of length i(ith order)
theta_mean = zeros(P,1);
% This is C_theta; identity vector of size i
theta_variance = 10000*eye(P);
% Posterior parameters
[big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
    theta_mean,theta_variance,y);

% Regenerate the audio
y_ML = G*theta_ML;
y_MAP = G*theta_MAP;
x = linspace(1, length(y), length(y));
x_fre = linspace(1, length(y)/2, length(y)/2);


% Plot power spectrum
power_origin = log(abs(fft(y, Fs1)).^2);
power_ML = log(abs(fft(y_ML, Fs1)).^2);
power_MAP = log(abs(fft(y_MAP, Fs1)).^2)
figure
plot(power_origin);
hold on
plot(power_ML);
hold on
plot(power_MAP);
xlim([8000 16000])
xlabel('Frequency')
ylabel('Amplitude')
legend('Original', 'ML', 'MAP')
title('Power Spectrum of AR Processes')
hold off