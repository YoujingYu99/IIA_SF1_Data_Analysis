clear all;
%add audio path
addpath('Week3_audio');
% Load both the defected and the original signals
[y_voice, Fs] = audioread('grosse_40_percent_missing.wav');
[y_orig, Fs] = audioread('grosse_original.wav');


y_original = y_voice;
% Initialise the same vector for output
y_out_ML = y_voice;
y_out_MAP = y_voice;
% Initialise model parameters
error_variance = 0.00000001;
resi_variance = 0;
alpha = 0.4;
N = length(y_original);

empty_signal_list = [];
% Check how many audio signals are lost
for i=1:N
    if y_original(i) == 0 && y_original(i+1) == 0 
        empty_signal_list = [empty_signal_list, i];
    end
end

% Get the starting and ending points of the lost packets
section_start_list = [empty_signal_list(1)];
section_end_list = [];
for j=1:length(empty_signal_list)-1
    % If a new block of empty signals
    if empty_signal_list(j+1) - empty_signal_list(j) > 1
        section_start_list = [section_start_list, empty_signal_list(j+1)]
        section_end_list = [section_end_list, empty_signal_list(j)]
    end
end

section_end_list = [section_end_list, empty_signal_list(end)]
        
% Generate different AR models for different sections
for j = 1:length(section_start_list)-1
    model_LLR_list = [];
    % Section before the pause
    % If first audio, se the start of music as position 1
    if j-1 <1
        before_block_start = 1;
    else
        before_block_start = section_end_list(j-1)+1;
    end
    % The section of audio before the silent block
    y_section = y_original(before_block_start :section_start_list(j)-1);
    N = length(y_section);
    % The section of audio after the silent block
    y_section_next = y_original(section_end_list(j)+1:section_start_list(j+1));
    length_of_empty_section = section_end_list(j) - section_start_list(j);
    % Find optimum model order
    for i = 1:20
        % Set different AR orders
        G_AR = ARmodel(y_section, N,i);
        G = G_AR;
        y = y_section(i+1:end);
        % ML estimation
        theta_ML = inv(transpose(G)*G)*transpose(G)*y;

        %MAP estimation
        % This is m_theta; column vector of length i(ith order)
        theta_mean = zeros(i,1);
        % This is C_theta; identity vector of size i
        theta_variance = eye(i);
        % Posterior parameters
        [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
            theta_mean,theta_variance,y);

        % Model Selection
        marginal_LLR = model_LLR(N,i,G,error_variance, theta_mean,theta_variance,y);
        model_LLR_list = [model_LLR_list; marginal_LLR];

    end

    % Return the optimum offset and the LLR
    [optimum_LLR, optimum_order] = max(model_LLR_list);
    % Convert back from the log domain to normal and normalise
    posterior_probs = exp(model_LLR_list);
    posterior_probs_norm = posterior_probs /sum(posterior_probs);

%     % Plot posterior probabilities against positions
%     figure(1)
%     n = linspace(1, 20, 20);
%     plot(n, posterior_probs_norm)
%     xlabel('theta')
%     ylabel('Probability Density')
%     title('Posterior Probabilities of Order of AR Process')

    P = optimum_order;
    % Regenerate the AR model
    G_AR = ARmodel(y_section, N,P);
    G = G_AR;
    y = y_section(P+1:end);
    % ML estimation
    theta_ML = inv(transpose(G)*G)*transpose(G)*y;

    %MAP estimation
    % This is m_theta; column vector of length i(ith order)
    theta_mean = zeros(P,1);
    % This is C_theta; identity vector of size i
    theta_variance = eye(P);
    % Posterior parameters
    [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
        theta_mean,theta_variance,y);

    % Forward and backward prediction
    y_missing_section_ML = zeros(length_of_empty_section,1);
    y_missing_section_forward_ML = zeros(length_of_empty_section,1);
    y_missing_section_backward_ML = zeros(length_of_empty_section,1);
    y_missing_section_MAP = zeros(length_of_empty_section,1);
    y_missing_section_forward_MAP = zeros(length_of_empty_section,1);
    y_missing_section_backward_MAP = zeros(length_of_empty_section,1);
    for k=1:length_of_empty_section 
        % Forward prediction
        y_missing_section_forward_ML(k) = fliplr(y_original(section_end_list(j)-P+k: section_end_list(j)+k-1))'*theta_ML + sqrt(resi_variance)*randn(1,1);
        y_missing_section_forward_MAP(k) = fliplr(y_original(section_end_list(j)-P+k: section_end_list(j)+k-1))'*theta_MAP + sqrt(resi_variance)*randn(1,1);
        % Backward prediction
        y_missing_section_backward_ML(k) = y_original(section_end_list(j)+k: section_end_list(j)+P+k-1)'*theta_ML + sqrt(resi_variance)*randn(1,1);
        y_missing_section_backward_MAP(k) = y_original(section_end_list(j)+k: section_end_list(j)+P+k-1)'*theta_MAP + sqrt(resi_variance)*randn(1,1);
        y_missing_section_ML(k) = alpha * y_missing_section_forward_ML(k) + (1-alpha)*y_missing_section_backward_ML(k);
        y_missing_section_MAP(k) = alpha * y_missing_section_forward_MAP(k) + (1-alpha)*y_missing_section_backward_MAP(k);

    end
    %Add the regenerated section back to original signal
    start = section_start_list(j):section_end_list(j);
    y_out_ML(section_start_list(j):section_end_list(j)-1) = y_missing_section_ML(:);
    y_out_MAP(section_start_list(j):section_end_list(j)-1) = y_missing_section_MAP(:);
end

% mse_ML = immse(y_orig, y_out_ML);
% mse_MAP = immse(y_orig, y_out_MAP);
% Plot the audio
figure(2)
plot(y_orig);
hold on
plot(y_out_ML);
hold on
plot(y_out_MAP);
xlabel('Time')
ylabel('Amplitude')
legend('Original', 'ML', 'MAP')
title('Waveform of Original and Estimated Processes')
hold off



% sound(x_noisy, Fs);
sound(y_voice, Fs);
% 
% % Sound the output
% sound(y_out_ML, Fs);
% rootdirectory = "C:\Users\Youjing Yu\Desktop\SF1\Inference\Processed_audio";
% file_name = fullfile(rootdirectory,'grosse_40_ML.wav');
% 
% audiowrite(file_name,y_out_ML', Fs);





