function [big_theta, big_phi, theta_MAP, theta_MAP_covariance] = posterior_terms(G,error_variance,...
    theta_mean,theta_variance,y)
% Compute big theta
big_theta = transpose(G)*y + error_variance*inv(theta_variance)*theta_mean;
% Compute big phi
big_phi = transpose(G)*G + error_variance*inv(theta_variance);
% Compute theta_MAP
theta_MAP = inv(big_phi)*big_theta;
theta_MAP_covariance = error_variance*inv(big_phi);
end
