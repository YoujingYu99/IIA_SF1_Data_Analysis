%add audio path
addpath('Processed_audio');

[x, Fs] = audioread('armst_37_missing_ML.wav');

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
noise_threshold = 2;
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
%listen to denoised signal
sound(y_out', Fs);

rootdirectory = 'C:\Users\Youjing Yu\Desktop'
file_name = fullfile( rootdirectory,'armst_37_missing_ML.wav');

audiowrite(file_name,y_out', Fs);
