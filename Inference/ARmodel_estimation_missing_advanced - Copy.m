clear;
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
resi_variance = 0.00000001;
alpha = 0.4;
N = length(y_original);

empty_signal_list = [];
% Check how many audio signals are lost
for i=1:N
    if y_original(i) == 0 & y_original(i+1) == 0 ;
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

section_end_list = [section_end_list, empty_signal_list(end)];

% Generate audios with the before, missing and after blocks.
% Create empty cell array
C = {length(section_start_list), 1};
% Need four points to decide: previous end, start, end, next start
% Special case for first section
indices_list = [1, section_start_list(1), section_end_list(1), section_start_list(2)];
C{1} = indices_list;
% Special case for last section
indices_list = [section_end_list(length(section_start_list) -1), section_start_list(length(section_start_list)),
    section_end_list(length(section_start_list)), length(y_original)];
C{end} = indices_list;
% For sections in the middle
for i = 2:length(section_start_list)-1
    indices_list = [section_end_list(i-1), section_start_list(i), section_end_list(i), section_start_list(i+1)];
    C{i} = indices_list;
end

% Generate different AR models for different sections
for k = 1:length(C)
    indices_list = C{k};
    x_a = y_original(indices_list(1):indices_list(2));
    xi = y_original(indices_list(2):indices_list(3));
    x_b = y_original(indices_list(3):indices_list(4));

    x_i = [x_a', x_b'];


end

% For all samples
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

    % Generate x_i
    indices_list = C{j};
    x_a = y_original(indices_list(1):indices_list(2));
    xi = y_original(indices_list(2)+1:indices_list(3));
    x_b = y_original(indices_list(3)+1:indices_list(4));

    x_i = [x_a', x_b'];

    % Generate the A matrix
    % The section length now depends on each section.
    section_length = indices_list(4) - indices_list(1)+1;
    A_ML = zeros(section_length-P, section_length);
    a_ML = [1, -1*theta_ML'];
    % Populate A_ML with estimated numbers
    for s=1:section_length-P
        A_ML(s,section_length-P+1-s:section_length-P-s+length(a_ML)) = a_ML;
    end
    A_ML = fliplr(A_ML);

    % Partition A
    A_ML_a = A_ML(:,1:length(x_a));
    Ai = A_ML(:,length(x_a)+1: length(x_a)+length(xi));
    A_ML_b = A_ML(:,length(x_a)+length(xi)+1: section_length);
    A_ML_i = [A_ML_a, A_ML_b]

    y_ML = -inv(Ai' * Ai)*Ai'* A_ML_i * x_i';

    % Generate the A matrix
    % The section length now depends on each section.
    A_MAP = zeros(section_length-P, section_length);
    a_MAP = [1, -1*theta_MAP'];
    % Populate A_MAP with estimated numbers
    for z=1:section_length-P
        A_MAP(z,section_length-P+1-z:section_length-P-z+length(a_MAP)) = a_MAP;
    end
    A_MAP = fliplr(A_MAP);

    % Partition A
    A_MAP_a = A_MAP(:,1:length(x_a));
    Ai = A_MAP(:,length(x_a)+1: length(x_a)+length(xi));
    A_MAP_b = A_MAP(:,length(x_a)+length(xi)+1: section_length);
    A_MAP_i = [A_MAP_a, A_MAP_b];


    y_MAP = -inv(Ai' * Ai)*Ai'* A_MAP_i * x_i';

    %Add the regenerated section back to original signal
    start = section_start_list(j):section_end_list(j);
    y_out_ML(section_start_list(j):section_end_list(j)-1) = y_ML(:);
    y_out_MAP(section_start_list(j):section_end_list(j)-1) = y_MAP(:);
end

sound(y_out_MAP, Fs);
rootdirectory = "C:\Users\Youjing Yu\Desktop\SF1\Inference\Processed_audio";
file_name1 = fullfile(rootdirectory,'grosse_40_advanced_MAP.wav');
file_name2 = fullfile(rootdirectory,'grosse_40_advanced_ML.wav');

audiowrite(file_name1,y_out_MAP', Fs);
audiowrite(file_name2,y_out_ML', Fs);
