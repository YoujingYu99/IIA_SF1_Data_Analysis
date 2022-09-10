clear all;
% Generate second order data
% Number of data points
N=100;
% Order of AR
P = 2;
% Set AR process parameters
error_variance = 10;
theta1 = 0.4;
theta2 = 0.35;
AR_model = arima('Constant',0,'AR',{theta1 theta2},'Variance',error_variance);
y_original = simulate(AR_model, 50);

% Generate second order data
% Number of data points
N=length(y_original);
order_list = [];

for j = 1:100
    variance_mag = 0.05*j;
    model_LLR_list = [];
    % MSE for different model orders
    for i = 1:20
        % Set different AR orders
        G_AR = ARmodel(y_original, N,i);
        G = G_AR;
        y = y_original(i+1:end);
        % ML estimation
        theta_ML = inv(transpose(G)*G)*transpose(G)*y;

        %MAP estimation
        % This is m_theta; column vector of length i(ith order)
        theta_mean = zeros(i,1);
        % This is C_theta; scaled identity vector of size i
        theta_variance = variance_mag * eye(i);
        % Posterior parameters
        [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
            theta_mean,theta_variance,y);

        % Model Selection
        marginal_LLR = model_LLR(N,i,G,error_variance, theta_mean,theta_variance,y);
        model_LLR_list = [model_LLR_list; marginal_LLR];

    end
  
    % Return the optimum offset and the LLR
    [optimum_LLR, optimum_order] = max(model_LLR_list);
    order_list = [order_list, optimum_order];

end


% Plot Model Order against prior variances
figure
n = linspace(1, 100, 100);
plot(n, order_list);
ylim([0,21]);
xlabel('Prior Variance')
ylabel('Optimum Order Estimated')
title('Optimum Order of AR Process against Prior Variances')





