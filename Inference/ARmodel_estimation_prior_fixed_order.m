clear all;
% Generate second order data
% Number of data points
N=100;
% Order of AR
P = 2;
% Set AR process parameters
error_variance = 1;
theta1 = 0.4;
theta2 = 0.35;
AR_model = arima('Constant',0,'AR',{theta1 theta2},'Variance',error_variance);
y_original = simulate(AR_model, 50);

% Generate second order data
% Number of data points
N=length(y_original);
ML_theta1_list = [];
ML_theta2_list = [];
MAP_theta1_list = [];
MAP_theta2_list = [];

for i = 1:100
    variance_mag = i;
    % Make AR model
    G_AR = ARmodel(y_original, N,2);
    G = G_AR;
    y = y_original(2+1:end);
    % ML estimation
    theta_ML = inv(transpose(G)*G)*transpose(G)*y;
    y_ML = G*theta_ML;
    ML_theta1_list = [ML_theta1_list, y_ML(1)];
    ML_theta2_list = [ML_theta2_list, y_ML(2)];

    %MAP estimation
    % This is m_theta; column vector of length i(ith order)
    theta_mean = zeros(2,1);
    % This is C_theta; identity vector of size i
    theta_variance = variance_mag * eye(2);
    % Posterior parameters
    [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
        theta_mean,theta_variance,y);
    y_MAP = G*theta_MAP;
    MAP_theta1_list = [MAP_theta1_list, y_MAP(1)];
    MAP_theta2_list = [MAP_theta2_list, y_MAP(2)];

end

figure(1)
plot(ML_theta1_list);
hold on
plot(MAP_theta1_list);
xlabel('Prior Variance')
ylabel('theta 1')
legend('ML', 'MAP')
title('Theta 1 estimation of ML and MAP with Increasing Prior Variances')
hold off

figure(2)
plot(ML_theta2_list);
hold on
plot(MAP_theta2_list);
xlabel('Prior Variance')
ylabel('theta 2')
legend('ML', 'MAP')
title('Theta 2 estimation of ML and MAP with Increasing Prior Variances')
hold off


figure(3)
plot(y);
hold on
plot(y_ML);
hold on
plot(y_MAP);
xlabel('Prior Variance')
ylabel('y Values')
legend('Original','ML', 'MAP')
title('Original and Estimated Processes using ML and MAP')
hold off





