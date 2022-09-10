clear all;

N=100;
%Order 1 here
P=1;
% Noise variance
sigma_n=2;
error_variance = sigma_n^2;
% Random signal
signal=[-3 5 -2 4 1 3 5 -1 2 4 6 5 -2 -2 1];


% Set prior mean to be 0
theta_mean = 0;

% Load data
str=load('hidden_data.mat');
y=str.y;

all_null_list = [];
offset_list = [];
for j=1:100
  % Prior variance
  sigma_theta=0.1*j;
  theta_variance = sigma_theta^2;
  % Create empty list for storing marginal likelihoods
  model_LLR_list=[];
  null_list = [];
  interval = [-0.001, 0.001];

  for i = 1:90
      % Generate empty array of length N
      G = zeros(N,1);
      % Add signal of length into array at each position, where i is the offset
      G(i:i+14) = signal';
      G = G(1:N,:);
      marginal_LLR = model_LLR(N,P,G,error_variance, theta_mean,theta_variance,y);
      model_LLR_list = [model_LLR_list; marginal_LLR];

      % Calculate theta_MAP
      [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
      theta_mean,theta_variance,y)

      % Calculate the probability that theta=0 using cdf
      null_prob = normcdf(interval, theta_MAP, sqrt(theta_variance));
      null_list = [null_list, null_prob(2)-null_prob(1)];
  end


  % Return the optimum offset and the LLR
  [optimum_LLR, optimum_offset] = max(model_LLR_list);
  offset_list = [offset_list, optimum_offset];
  % Convert back from the log domain to normal and normalise
  posterior_probs = exp(model_LLR_list);
  posterior_probs_norm = posterior_probs /sum(posterior_probs);



  % Now that we know the position of offset, recalculate ML and MAP solution
  % Regenerate empty array of length N
  G = zeros(N,1);
  % Recover the signal
  G(optimum_offset:optimum_offset+14) = signal';
  % Calculate theta_ML and its variance
  theta_ML = inv(transpose(G)*G)*transpose(G)*y;
  theta_ML_variance = error_variance*inv((transpose(G)*G));
  % Calculate theta MAP using the optimal G
  [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
      theta_mean,theta_variance,y);

  % Plot ML and MAP probabilities
  x = [-1.5:0.001:1.5];
  y_prior = normpdf(x,theta_mean, sqrt(theta_variance));
  y_likelihood = normpdf(x,theta_ML,sqrt(theta_ML_variance));
  y_posterior = normpdf(x,theta_MAP,sqrt(theta_MAP_covariance));

  % Calculate the null probability
  null_prob_final = mean(null_list);
  all_null_list = [all_null_list, null_prob_final];
end



figure(1)
n = linspace(0.1, 10, 100);
plot(n.^2, all_null_list);
xlabel('Prior Variance')
ylabel('Null Probability')
title('Null Probability against Prior Variance')

figure(2)
n = linspace(0.1, 10, 100);
plot(n.^2, offset_list);
xlabel('Prior Variance')
ylabel('Optimual Offset')
title('Optimual Offset against Prior Variance')
