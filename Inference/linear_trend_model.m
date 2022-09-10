% Set theta and N
% True mean
theta1=5;
theta2=10;
theta = [theta1; theta2];
% theta = theta';
% Number of data points
N=100;
% Variance of error
error_variance = 1000;

%Create N*1 matrix of ones for g1
g1 = ones(N,1);
% Create matrix of 1 to N for g2
g2 = [1:N]';
% Concatenate; P=2.
G = [g1 g2];
y = G*theta + sqrt(error_variance)*randn(N,1);
y_theory = G*theta;
figure(1)
plot(y);
hold on
plot(y_theory);
xlabel('n')
ylabel('y_n')
legend('Generated Data','Theoretical Distribution')
title(['Data with theta1=5, theta2=10, error variance=1000, N=', num2str(N)])

% ML estimate for theta
theta_ML = inv(transpose(G)*G)*transpose(G)*y;

% Calculated from equation (19) when C_theta tends to infinity
theta_ML_variance = error_variance*inv((transpose(G)*G));



% MAP estimate for theta

% Set theta prior parameters
% This is m_theta
theta_mean = [5;10];
% This is C_theta
theta_variance = 10*eye(2);


% Posterior parameters
[big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
    theta_mean,theta_variance,y);

% Plot ML and MAP methods for theta1
x = [3:0.001:7];
y_prior_1 = normpdf(x,theta_mean(1), sqrt(theta_variance(1,1)));
y_likelihood_1 = normpdf(x,theta_ML(1),sqrt(theta_ML_variance(1,1)));
y_posterior_1 = normpdf(x,theta_MAP(1),sqrt(theta_MAP_covariance(1,1)));
figure(2)
plot(x, y_prior_1)
plot(x, y_likelihood_1)
plot(x, y_posterior_1)
xline(theta_ML(1),'r','ML','LabelHorizontalAlignment','left')
xline(theta_MAP(1),'m','MAP','LabelHorizontalAlignment','right')
legend('Prior','Likelihood','Posterior')
xlabel('theta')
ylabel('Probability Density')
title(sprintf('The Probability Distributions with Prior 1 Mean=1 and Prior Variance = 1'))

% Plot ML and MAP methods for theta2
x = [9.8:0.001:10.2];
y_prior_2 = normpdf(x,theta_mean(2), sqrt(theta_variance(2,2)));
y_likelihood_2 = normpdf(x,theta_ML(2),sqrt(theta_ML_variance(2,2)));
y_posterior_2 = normpdf(x,theta_MAP(2),sqrt(theta_MAP_covariance(2,2)));
figure(3)
plot(x, y_prior_2)
plot(x, y_likelihood_2)
plot(x, y_posterior_2)
xline(theta_ML(2),'r','ML','LabelHorizontalAlignment','left')
xline(theta_MAP(2),'m','MAP','LabelHorizontalAlignment','right')
legend('Prior','Likelihood','Posterior')
xlabel('theta')
ylabel('Probability Density')
title(sprintf('The Probability Distributions with Prior 2 Mean=2 and Prior Variance = 1'))
