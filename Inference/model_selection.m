clear all;
% First model
theta1=0;
% Second model
theta2 = 10;
% Third model
theta3_1= 10;
theta3_2 = 1;
theta3 = [theta3_1; theta3_2];
theta_variance_prior = 1;
% Number of data points
N=100;
P = 0;
% Variance of error
error_variance = 100000;

% Generate datasets
% Pure noise
G1 = 0;
y1 = sqrt(error_variance)*randn(N,1);

% Linear model
G2 = ones(N,1);
% Generate data
y2 = G2*theta2 + sqrt(error_variance)*randn(N,1);


% Linear trend
%Create N*1 matrix of ones for g1
g1 = ones(N,1);
% Create matrix of 1 to N for g2
g2 = [1:N]';
% Concatenate; P=2.
G3 = [g1 g2];
y3 = G3*theta3 + sqrt(error_variance)*randn(N,1);


theta_mean_1 = 0;
theta_mean_2 = 0;
theta_mean_3 = [0;0];


% % theta 1
% theta1_array_ML = zeros(1,10);
% theta1_array_MAP = zeros(1,10);
% % theta 2
% theta2_array_ML = zeros(1,10);
% theta2_array_MAP = zeros(1,10);
% 
% % Try for various prior variances
% for i=1:10
% MAP array
MAP_array_1 = zeros(1,10);
MAP_array_2 = zeros(1,10);
MAP_array_3 = zeros(1,10);

% Try for various prior variances
for i=1:10000
    theta_variance_prior = 10*i;
    theta_variance_2 = theta_variance_prior;
    theta_variance_3 = theta_variance_prior*eye(2);
    

    % MAP estimation
    M1_marginal_1 = model_LLR(N,P,G1,error_variance, theta_mean_1,theta_variance_prior,y3);
    M1_marginal_2 = model_LLR(N,P,G2,error_variance, theta_mean_2,theta_variance_2,y3);
    M1_marginal_3 = model_LLR(N,P,G3,error_variance, theta_mean_3,theta_variance_3,y3);
    
    %     Record result in arrays
    MAP_array_1(i) = M1_marginal_1;
    MAP_array_2(i) = M1_marginal_2;
    MAP_array_3(i) = M1_marginal_3;
end




figure(3)
plot(MAP_array_1);
hold on
plot(MAP_array_2);
hold on
plot(MAP_array_3);
% hold on
xlabel('Prior Variance')
ylabel('Marginal Probabilities')
legend('Pure Noise Model', 'Linear Model', 'Linear Trend Model')
title('Probabilities of Models on Data Generated from Model 3')
hold off




