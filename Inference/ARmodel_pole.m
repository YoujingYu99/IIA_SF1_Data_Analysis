clear all;

% Generate second order data
% Number of data points
N=100;
% Order of AR
P = 2;
% Set AR process parameters
error_variance = 1;
theta1 = 0.4;
theta2 = 1.2;

% 3F3 Notes
% Number of data points
N=1000;
% AR order
P=2;
pole(1)=theta2*exp(j*0.1*pi);
pole(2)=conj(pole(1));
% Coefficients of a polynomial whose roots are the poles
a=poly(pole);
sigma_e=1;
e=randn(N,1)*sigma_e;
% Past data points assumed to be zero:
x=filter(1,a,e);
figure(1)
subplot(211), plot(e)
title('e_n')
subplot(212), plot(x)
title('x_n')

figure(10);
plot(x)
title('Simulated AR(2) Process with a1=0.4, a2=1.2')
% Make a Matlab system model:
sys=tf(1,a,1);
figure, pzplot(sys)
title('Poles of AR(2) Model with a1=0.4, a2=1.2')
[H,w]=freqz(1,a);
figure(3)
subplot(211),
semilogy(w,abs(H).^2*sigma_e^2)
title('Power Spectrum of AR(2) Model with a1=0.4, a2=1.2')
subplot(212),
X=(abs(fft(x)).^2);
semilogy((0:N/2-1)*pi*2/N,X(1:N/2))
title('|DFT|^2 of x_n')


AR_model = arima('Constant',0,'AR',{theta1 theta2},'Variance',error_variance);
y = simulate(AR_model, 50);
figure(4)
plot(y)
xlim([0,50])
title('Simulated AR(2) Process with a1=0.4, a2=1.2')



