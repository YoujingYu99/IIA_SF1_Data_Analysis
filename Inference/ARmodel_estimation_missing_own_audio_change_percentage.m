clear all;
%add audio path
addpath('Week3_audio');
% Load both the defected and the original signals

[y_orig, Fs] = audioread('armst_37_orig.wav');

% Vary percentage of data missing; section_length constant at 1000
ML_MSE_list1 = [];
MAP_MSE_list1 = [];
for i=1:20
    percentage = i*0.01;
    [mse_ML, mse_MAP] = AR_mse(y_orig, Fs, percentage, 1000, 0);
    ML_MSE_list1 = [ML_MSE_list1, mse_ML];
    MAP_MSE_list1 = [MAP_MSE_list1, mse_MAP];
end

% Plot MSE when varing percentage
figure(1)
n = linspace(1, 20, 20);
plot(n, ML_MSE_list1);
hold on
plot(n, MAP_MSE_list1);
xlabel('Percentage of Data Missing')
ylabel('Mean Squared Error')
legend('ML', 'MAP')
title('MSE against % of Missing Data with ML and MAP')



% 
% % Vary length of section; percentage constant at 0.1
% ML_MSE_list2 = [];
% MAP_MSE_list2 = [];
% for i=1:200
%     section_length = 50*i;
%     [mse_ML, mse_MAP] = AR_mse(y_orig, Fs, 0.1, section_length, 0);
%     ML_MSE_list2 = [ML_MSE_list2, mse_ML];
%     MAP_MSE_list2 = [MAP_MSE_list2, mse_MAP];
% end
% 
% % Plot MSE when varing section_length
% figure(2)
% n = linspace(1, 200, 200);
% plot(n*50, ML_MSE_list2);
% hold on
% plot(n*50, MAP_MSE_list2);
% xlabel('Section Length')
% ylabel('Mean Squared Error')
% legend('ML', 'MAP')
% title('MSE against Section Length with ML and MAP')

% % Vary residual noise variance
% ML_MSE_list3 = [];
% MAP_MSE_list3 = [];
% for i=1:200
%     resi_variance = i*0.01;
%     [mse_ML, mse_MAP] = AR_mse(y_orig, Fs, 0.1, 1000, resi_variance);
%     ML_MSE_list3 = [ML_MSE_list3, mse_ML];
%     MAP_MSE_list3 = [MAP_MSE_list3, mse_MAP];
% end
% 
% % Plot MSE when varing percentage
% figure(1)
% n = linspace(1, 200, 200);
% plot(n*0.01, ML_MSE_list3);
% hold on
% plot(n*0.01, MAP_MSE_list3);
% xlabel('Residual Noise Variance')
% ylabel('Mean Squared Error')
% legend('ML', 'MAP')
% title('MSE for Residual Noise Variance with ML and MAP')







