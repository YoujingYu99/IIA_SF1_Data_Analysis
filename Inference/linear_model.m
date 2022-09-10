% Set theta and N
% True mean
theta=5;
% Number of data points
N=50;
% Variance of error
error_variance = 10;

% Create N*1 matrix for G
G = ones(N,1);
% Generate data
y = G*theta + sqrt(error_variance)*randn(N,1);
figure(1)
plot(y)
yline(theta, 'r', 'theta')
xlabel('n')
ylabel('y_n')
legend('Generated Data', 'Theoretical Distribution')
title(['Example data with theta=5 and N=', num2str(N)])

% ML estimate for theta
% theta_ML = inv(transpose(G)*G)*transpose(G)*y;
theta_ML = mean(y);
% Calculated from equation (19) when C_theta tends to infinity
theta_ML_variance = error_variance/(transpose(G)*G);



% MAP estimate for theta

% Set theta prior parameters
% This is m_theta
theta_mean = 5;
% This is C_theta
theta_variance = 10;


% Posterior parameters
[big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
    theta_mean,theta_variance,y);

% Plot ML and MAP methods
x = [3:0.001:7];
y_prior = normpdf(x,theta_mean, sqrt(theta_variance));
y_likelihood = normpdf(x,theta_ML,sqrt(theta_ML_variance));
y_posterior = normpdf(x,theta_MAP,sqrt(theta_MAP_covariance));
figure(2)
plot(x, y_prior)
plot(x, y_likelihood)
plot(x, y_posterior)
xline(theta_ML,'r','ML','LabelHorizontalAlignment','left')
xline(theta_MAP,'m','MAP','LabelHorizontalAlignment','right')
legend('Prior','Likelihood','Posterior')
xlabel('theta')
ylabel('Probability Density')
title(sprintf('The Parameter Distributions with Prior Variance=10'))
