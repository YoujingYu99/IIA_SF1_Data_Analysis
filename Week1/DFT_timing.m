%set length to be investigated
length_array = linspace(256, 256*4, 256*4+1);
length_of_array = length(length_array);
%choose a given audio to analyse
signal = audioread('armst_37_orig.wav');
%initialise timings array for fft and dft
times_fft=zeros(1,length(length_array));
times_dft=zeros(1,length(length_array));

%time dft calculation
tic;
for i = 1:length_of_array
    N_chosen = length_array(i);
    %calculate DFT
    signal_dft = DFT_vectorized(signal, N_chosen);
    times_dft(i) = toc;
end
%taking natural log of the timings data
log_times_dft = log(times_dft);

%time fft calculation
tic;
for i = 1:length_of_array
    N_chosen = length_array(i);
    %calculate FFT
    signal_fft = fft(signal, N_chosen);
    times_fft(i) = toc;
end
%taking natural log of the timings data
log_times_fft = log(times_fft);

% %plot timings in log domain
% subplot(1,2,1);
% plot(length_array, log_times_dft);
% title('Timing for DFT Function')
% xlabel('N')
% ylabel('Time(sec)(log)')
% subplot(1,2,2);
% plot(length_array, log_times_fft);
% title('Timing for FFT Function')
% xlabel('N')
% ylabel('Time(sec)(log)')

%plot complexity in log domain
subplot(1,2,1);
plot(length_array, log(length_array.^2));
title('Complexity for DFT Function')
xlabel('N')
ylabel('Relative Time(log)')
subplot(1,2,2);
plot(length_array, log(log(length_array).* length_array));
title('Complexity for FFT Function')
xlabel('N')
ylabel('Relative Time(log)')
