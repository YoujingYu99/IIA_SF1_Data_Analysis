function  model_selection_plot(original_data_model) 
        % Set parameters
        % First model
        theta1=0;
        % Second model
        theta2 = 5;
        % Third model
        theta3_1=10;
        theta3_2 = 0.1;
        theta3 = [theta3_1; theta3_2];
        theta_variance_prior = 100;
        % Number of data points
        N=100;
        % Variance of error
        error_variance = 1000;
        
        G1 = 0;
        G2 = ones(N,1);
        % Linear trend
        %Create N*1 matrix of ones for g1
        g1 = ones(N,1);
        % Create matrix of 1 to N for g2
        g2 = [1:N]';
        % Concatenate; P=2.
        G3 = [g1 g2];
    % Generate y values
    if original_data_model == 1
        y = sqrt(error_variance)*randn(N,1);
    elseif original_data_model == 2
        % Generate data
        y = G2*theta2 + sqrt(error_variance)*randn(N,1);
    else
        y = G3*theta3 + sqrt(error_variance)*randn(N,1);
    end


    % Set theta prior parameters
    % This is m_theta_1
    theta_mean_1 = 0;
    % This is C_theta_1
    theta_variance_1 = 0;

    % Posterior parameters
    [big_theta_1, big_phi_1, theta_MAP_1, theta_MAP_covariance_1] = posterior_terms(G1,error_variance,...
        theta_mean_1,theta_variance_1,y);

    % Set theta prior parameters
    % This is m_theta_2
    theta_mean_2 = 0;
    % This is C_theta_2
    theta_variance_2 = theta_variance_prior;

    % Posterior parameters
    [big_theta_2, big_phi_2, theta_MAP_2, theta_MAP_covariance_2] = posterior_terms(G2,error_variance,...
        theta_mean_2,theta_variance_2,y);

    % Set theta prior parameters
    % This is m_theta_3
    theta_mean_3 = [0;0];
    % This is C_theta_3
    theta_variance_3 = theta_variance_prior*eye(2);

    % Posterior parameters
    [big_theta_3, big_phi_3, theta_MAP_3, theta_MAP_covariance_3] = posterior_terms(G3,error_variance,...
        theta_mean_3,theta_variance_3,y);

    % Regenerate the data from 3 models
    y_model1 = G1*theta_MAP_1;
    y_model2 = G2*theta_MAP_2;
    y_model3 = G3*theta_MAP_3;

    figure(1)
    plot(y);
    hold on
    plot(y_model1);
    hold on
    plot(y_model2);
    hold on
    plot(y_model3);
    xlabel('Sequence')
    ylabel('Data')
    legend('Original','Constant Model', 'Linear Model', 'Linear Trend Model')
    title('Original and Data Regenerated from Models 1, 2,3')
    hold off
end
