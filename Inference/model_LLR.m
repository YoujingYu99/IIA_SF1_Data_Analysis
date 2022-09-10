function marginal_LLR = model_LLR(N,P,G,error_variance, theta_mean,theta_variance,y)
% Compute the posterior terms first
[big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
    theta_mean,theta_variance,y);
exp_term = -(transpose(y)*y + error_variance*transpose(theta_mean)*inv(theta_variance)*theta_mean - transpose(big_theta)*theta_MAP)/(2*error_variance);
norm_term = -P*log(2*pi)/2 - log(det(theta_variance)*det(big_phi))/2 -(N-P)*log(2*pi*error_variance)/2;
marginal_LLR = norm_term + exp_term;
disp(norm_term);
end
