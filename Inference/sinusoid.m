clear all;

N = 50;
% Order of model
P = 1;
G = zeros(N, 2*P);
% Set bins of frequencies
w_bin = linspace(0, pi, 1);


% Populate the model
for i =1:N
    for j = 1:P
        % Assign the cosines and sines
        G(i, 2*j-1) = cos(i*w_bin(j));
        G(i, 2*j) = sin(i*w_bin(j));
    end
end

% Data generation
a1 = 0.3;
b1 = 0.4;
theta = [a1, b1]';
error_variance = 0.001;
y = G*theta + sqrt(error_variance)*randn(N,1);


% ML estimation
theta_ML = inv(transpose(G)*G)*transpose(G)*y;
% MAP estimation
% Set theta prior parameters
% This is m_theta
theta_mean = [0; 0];
% This is C_theta
theta_variance = 10*eye(2);
[big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
    theta_mean,theta_variance,y);

% Estimate the output
y_ML = G*theta_ML;
y_MAP = G*theta_MAP;

% Plot results to see
plot(y);
hold on
plot(y_ML);
hold on
plot(y_MAP);
hold on 
xlabel('Time')
ylabel('Data')
legend('Original', 'ML', 'MAP')
title('Waveform of Original and Estimated Processes using Sinusoidal Model Order=2')
hold off

fv_plot = fvtool(y, 1, y_ML, 2, y_MAP, 3);
legend(fv_plot, 'Original', 'ML', 'MAP');






