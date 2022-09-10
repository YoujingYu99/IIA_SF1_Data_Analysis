clear;

%add audio path
addpath('Audio_Week2');
[x, Fs] = audioread('grosse_original.wav');

%noise addition
noise_amp = 0.01;


%set length and mse array
length_array = linspace(256, 256*12, 12);
mse_array_spec = zeros(length(length_array), 1);
mse_array_wiener = zeros(length(length_array), 1);
mse_array_power = zeros(length(length_array), 1);

for i = 1:length(length_array)
    N=length_array(i);
    overlap=N/2;

    %calulate signal with noise, denoised signal and MSE with filter type chosen
    [x_noisy, y_out, mse1] = own_filter_func(N, overlap, x, Fs, noise_amp,'white' ,'spectralsub');
    mse_array_spec(i) = mse1;
    [x_noisy, y_out, mse2] = own_filter_func(N, overlap, x, Fs, noise_amp,'white', 'Wiener');
    mse_array_wiener(i) = mse2;
    [x_noisy, y_out, mse3] = own_filter_func(N, overlap, x, Fs, noise_amp, 'white','powersub');
    mse_array_power(i) = mse3;
end

N=1792;
overlap=N/2;
[x_noisy, y_out, mse] = own_filter_func(N, overlap, x, Fs, noise_amp, 'white', 'spectralsub');

sound(y_out', Fs);

figure;
set(gcf,'position',[0, 0, 1000, 500]);
subplot(1,2,1);
plot(x_noisy);
title('Time Response Before Filtering')
xlabel('Time')
ylabel('Magnitude')

subplot(1,2,2);
plot(y_out);
title('Time Response After Filtering')
xlabel('Time')
ylabel('Magnitude')

figure;
set(gcf,'position',[0, 0, 1000, 500]);
subplot(1,2,1);
plot(log(abs(fft(x_noisy))));
title('Frequency Response Before Filtering')
xlabel('Frequency')
ylabel('Magnitude(log)')

subplot(1,2,2);
plot(log(abs(fft(y_out))));
title('Frequency Response After Filtering')
xlabel('Frequency')
ylabel('Magnitude(log)')





