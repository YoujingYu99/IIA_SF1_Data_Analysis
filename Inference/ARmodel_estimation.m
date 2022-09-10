clear;
% Generate second order data
% Number of data points
N=100;
% Order of AR
P = 2;
% Set AR process parameters
error_variance = 100;
theta1 = 0.4;
theta2 = 0.35;
AR_model = arima('Constant',0,'AR',{theta1 theta2},'Variance',error_variance);
y_original = simulate(AR_model, 50);

% Generate second order data
% Number of data points
N=length(y_original);



ML_MSE_list = [];
MAP_MSE_list = [];
model_LLR_list = [];

% MSE for different model orders
for i = 1:20
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
    theta_variance = 0.1*eye(i);
    % Posterior parameters
    [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
        theta_mean,theta_variance,y);
    y_MAP = G*theta_MAP;
    MAP_MSE = mean((y_MAP - y).^2);
    MAP_MSE_list = [MAP_MSE_list, MAP_MSE];
    
    % Model Selection
    marginal_LLR = model_LLR(N,i,G,error_variance, theta_mean,theta_variance,y);
    model_LLR_list = [model_LLR_list; marginal_LLR];
    
end
  
figure(1)
plot(ML_MSE_list);
hold on
plot(MAP_MSE_list);
hold on
xlabel('Model Order')
ylabel('Mean Squared Error')
legend('ML', 'MAP')
title('MSE for ML and MAP Estimation Against Model Order with P=2')


% Return the optimum offset and the LLR
[optimum_LLR, optimum_order] = max(model_LLR_list);
% Convert back from the log domain to normal and normalise
posterior_probs = exp(model_LLR_list);
posterior_probs_norm = posterior_probs /sum(posterior_probs);


% Plot posterior probabilities against positions
figure(2)
n = linspace(1, 20, 20);
plot(n, posterior_probs_norm)
xlabel('Model Order')
ylabel('Probability Density')
title('Posterior Probabilities of Order of AR Process')



