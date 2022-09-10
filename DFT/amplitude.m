%set length to be investigated
length_array = linspace(256, 256*32, 256*31+1);
length_of_array = length(length_array);
normalised_freq = 0.5;
N = 256;
n = 1:N;
%create the signal required
signal = exp(j*normalised_freq*pi*n);


%different forms of signals
%linear increase
A = 0.5;
B = 5;
linear_gain = A*n + 1
sig1 = linear_gain.*signal

%periodic modulation
phi = 0.01;
beta = 0.05;
periodic_mod = 1 + beta*(sin(phi*n));
sig2 = periodic_mod.*signal;

%AR modulation
AR_mod = zeros(1,N);
AR_mod(1) = 1; 
for i = 2:N - 1
    v_n = 0.05*rand(1,1);
    AR_mod(i+1) = AR_mod(i) + v_n;
end
sig3 = AR_mod.*signal;
    


fv_am = fvtool(signal, 1,sig1, 2, sig2, 3, sig3, 4);
legend(fv_am, 'original', 'Linear Modulation', 'Periodic modulation', 'Autoregressive');
