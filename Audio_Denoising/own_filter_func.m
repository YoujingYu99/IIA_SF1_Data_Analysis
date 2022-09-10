function [x_noisy, y_out, mse] = own_filter_func(N, overlap, x, Fs, noise_amp, noise_type, noise_reduction)

    %get length of data
    signal_length = length(x);
    if strcmpi(noise_type, 'white')  
        %white noise addition
        n1 = noise_amp*randn(signal_length,1)';
    elseif strcmpi(noise_type, 'pink')       
        %pink noise addition
    %     n1 = noise_amp* pinknoise(signal_length, 'like',x);
        cn = dsp.ColoredNoise(1, signal_length ,1);
        rng default;
        n1 = noise_amp*cn()';
        disp(size(n1));
    end

    x_noisy = x' + n1;

    %calculate power spectrum of noise
    Sn = mean(abs(fft(n1, Fs)).^2);
    disp(Sn);

    %new array of same length as the signal, which is to be the output
    y_out=0*x_noisy;
    %partition x into data segments of length N and overlap length N/2
    x_new=buffer(x_noisy,N,overlap);

    %set number of samples and number of frames
    [N_samps, N_frames]=size(x_new);
    %create a matrix of N_frames columns.
    %each column is a hanning window of length N so the matrix is N*N_frames
    x_w=repmat(hanning(N),1,N_frames).*x_new;

    for frame_no=1:N_frames-2
        %fft on each column of hanning window
        X_w(:,frame_no)=fft(x_w(:,frame_no));

        %duplicate into a second matrix
        Y_w(:,frame_no)=X_w(:,frame_no);
        
        %Wiener filter
        z = Y_w(:, frame_no);
        SNR = (norm(z)^2 - Sn)/Sn;
        if strcmpi(noise_reduction, 'Wiener')            
            if SNR > 0
                f_z = SNR/(1+SNR); 
            else
                f_z = 0;
            Y_w(:, frame_no) = f_z*Y_w(:, frame_no);
            end
        elseif strcmpi(noise_reduction, 'spectralsub') 
            if SNR > 0
                f_z = 1 - 1/(sqrt(1+SNR)); 
            else
                f_z = 0;
            Y_w(:, frame_no) = f_z*Y_w(:, frame_no);
            end
        elseif strcmpi(noise_reduction, 'powersub') 
            if SNR > 0
                f_z = sqrt(SNR/(1+SNR)); 
            else
                f_z = 0;
            Y_w(:, frame_no) = f_z*Y_w(:, frame_no);
            end
        end
        %compute ifft of Y_w and copy into y_w
        y_w(:,frame_no)=ifft(Y_w(:,frame_no));

        %final result obtained by concatenating the frames and throwing away
        %the overlaps
        y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+y_w(:,frame_no)';
    end    

    %normalise the output between 0 and 1
    y_out_norm = y_out/max(y_out);
    x_norm = x'/max(x);
    
    %calculate mse
    mse = immse(x_norm, y_out_norm);
        