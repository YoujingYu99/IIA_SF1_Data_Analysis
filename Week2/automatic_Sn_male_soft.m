%add audio path
addpath('Noisy_audio');

[x, Fs] = audioread('male_speech_noisy_soft.wav');
% sound(x, Fs)

%noise psd estimation
%set data length N and blocklength which is half of the data length
N=2048;
overlap=N/2;

%calulate signal with noise, denoised signal and MSE with filter type chosen
x_noisy = x';
figure;
plot(abs(x_noisy));

%calculate power spectral density of noise automatically
x_clean = abs(x_noisy);
figure;
%plot the time response to determine a threshold 
plot(x_clean);
threshold = 0.3;
%set the first and last example to be 0
x_clean(1) = 0;
x_clean(end) = 0;
for i=2:length(x_noisy)-1
    %if the max in three samples is smaller than threshold, deem as noise.
    %otherwise keep.
    if max([x_clean(i-1), x_clean(i), x_clean(i+1)]) < threshold
        x_clean(i) = 0;
    end
end
%subtract the arrays to determine the noise array.
x_noise = abs(x_noisy) - x_clean;
x_noise_fft = abs(fft(x_noise));
Sn = mean(x_noise_fft.^2);
    



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

    %Wiener filter
    z = Y_w(:, frame_no);
    SNR = (norm(z)^2 - Sn)/Sn;
    % Wiener       
    if SNR > 0
        f_z = SNR/(1+SNR); 
    else
        f_z = 0;
    Y_w(:, frame_no) = f_z*Y_w(:, frame_no);
    end
%   'spectralsub'
%     if SNR > 0
%         f_z = 1 - 1/(sqrt(1+SNR)); 
%     else
%         f_z = 0;
%     Y_w(:, frame_no) = f_z*Y_w(:, frame_no);
%     end
%     elseif strcmpi(noise_reduction, 'powersub') 
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

% sound(x_noisy, Fs);
sound(y_out', Fs);
figure;
subplot(1,2,1);
plot(abs(fft(x)));
subplot(1,2,2);
plot(abs(fft(y_out)));



