[x, Fs] = audioread('grosse_original.wav');
N = length(x)

cn = dsp.ColoredNoise(1, Fs ,1);
rng default;
n_pink = cn();
fft_pink = log(abs(fft(n_pink)));

n_white = randn(N,1)';
fft_white = log(abs(fft(n_white)));



figure;
set(gcf,'position',[0, 0, 1000, 500]);
subplot(1,2,1);
plot(fft_pink);
title('DFT of Pink Noise')
xlabel('Frequency')
ylabel('Magnitude(log)')

subplot(1,2,2);
plot(fft_white);
title('DFT of White Noise')
xlabel('Frequency')
ylabel('Magnitude(log)')
