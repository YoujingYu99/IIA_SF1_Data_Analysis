signal = audioread('armst_37_orig.wav');
N_chosen = 1024;
signal_fft = fft(signal, N_chosen);
signal_dft = DFT_vectorized(signal, N_chosen);
subplot(1,2,1);
plot(abs(signal_fft));
subplot(1,2,2);
plot(abs(signal_dft));



