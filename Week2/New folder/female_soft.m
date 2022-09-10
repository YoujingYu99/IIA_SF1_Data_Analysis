%add audio path
addpath('Noisy_audio');

[x, Fs] = audioread('female_speech_noisy_soft.wav');

%set data length N and blocklength to be half of the data length
N=2048;
overlap=N/2;

%calulate signal with noise, denoised signal and MSE with filter type chosen
x_noisy = x';

%calculate power spectral density of noise
x_fft = abs(fft(x_noisy));
plot(x_fft);
disp(max(x_fft));
x_noise = x_fft;
noise_threshold = 1.0;
%discard signals. Only retain noise
x_noise(x_noise < noise_threshold) = [];
Sn = mean(x_noise.^2);

%new array of same length as the signal, which is to be the output
y_out=0*x_noisy;
%partition x into data segments of length N and overlap length N/2
x_new=buffer(x_noisy,N,overlap);

%set number of samples and number of frames
[N_samps,N_frames]=size(x_new);

%create a matrix of N_frames columns.
%each column is a hanning window of length N so the matrix is N*N_frames
x_w=repmat(hamming(N),1,N_frames).*x_new;

for frame_no=1:N_frames-2
    %fft on each column of hanning window
    X_w(:,frame_no)=fft(x_w(:,frame_no));

    %duplicate into a second matrix
    Y_w(:,frame_no)=X_w(:,frame_no);

    %Calculate SNR
    z = Y_w(:, frame_no);
    SNR = (norm(z)^2 - Sn)/Sn;
    % Wiener       
    if SNR > 0
        f_z = SNR/(1+SNR); 
    else
        f_z = 0;
    Y_w(:, frame_no) = f_z*Y_w(:, frame_no);
    end
%   %spectral subtraction
%     if SNR > 0
%         f_z = 1 - 1/(sqrt(1+SNR)); 
%     else
%         f_z = 0;
%     Y_w(:, frame_no) = f_z*Y_w(:, frame_no);
%     end
%   %power subtraction 
%         if SNR > 0
%             f_z = sqrt(SNR/(1+SNR)); 
%         else
%             f_z = 0;
%         Y_w(:, frame_no) = f_z*Y_w(:, frame_no);
%         end
    %compute ifft of Y_w and copy into y_w
    y_w(:,frame_no)=ifft(Y_w(:,frame_no));

    %final result obtained by concatenating the frames and throwing away
    %the overlaps
    y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+y_w(:,frame_no)';
end    

%normalise the output
y_out_norm = y_out/max(y_out);
x_norm = x'/max(x);

%calculate mse
mse = immse(x_norm, y_out_norm);

%listen to denoised signal
sound(y_out', Fs);

%plot time and frequency response
figure;
set(gcf,'position',[0, 0, 1000, 500]);
subplot(1,2,1);
plot(log(abs(fft(x))));
title('Frequency Response Before Filtering')
xlabel('Frequency')
ylabel('Magnitude(log)')

subplot(1,2,2);
plot(log(abs(fft(y_out))));
title('Frequency Response After Filtering')
xlabel('Frequency')
ylabel('Magnitude(log)')

figure;
set(gcf,'position',[0, 0, 1000, 500]);
subplot(1,2,1);
plot(x);
title('Time Response Before Filtering')
xlabel('Time')
ylabel('Magnitude')

subplot(1,2,2);
plot(y_out);
title('Time Response After Filtering')
xlabel('Time')
ylabel('Magnitude')


% rootdirectory = '\\cued-fs\users\General\yy471\windows-home\Desktop\SF1\Week2\New folder'
% file_name = fullfile( rootdirectory,'female_soft.wav');
% 
% audiowrite(file_name,y_out', Fs);
