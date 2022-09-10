N=100;
% Noise variance
sigma_n=2;
% Prior variance
sigma_theta=1;
% Random signal
signal=[-3 5 -2 4 1 3 5 -1 2 4 6 5 -2 -2 1];

% Set noise and prior parameters
noise=sigma_n*randn(N,1);
theta=sigma_theta*randn(1);
% 90 uniformly distributed random numbers which are rounded off
offset=round(90*rand(1));

y=noise;
signal_offset=0*noise;
% Choose the place to hide data
signal_offset(offset:offset+14)=signal’*theta;
y=y+signal_offset;
