function audio_window(signal, window_type, window_length, Fs)
%# Plot the single-sided frequency response of a signal convolved with a
%# chosen window in the log domain, with a specified window length and 
%# sampling frequency.
%specifying window types
if window_type == 'hamm'
    signal = signal.* hamming(window_length);
elseif window_type == 'hann'
    signal = signal.* hanning(window_length);
elseif window_type == 'bartlett'
    signal = signal.* bartlett(window_length);
elseif window_type == 'kaiser'
    signal = signal.* kaiser(window_length);
end

%plot single-sided spectrum
signal_fft = fft(signal);
L = length(signal_fft);

P2 = abs(signal_fft/L);
%only retaining one side of information
P1 = P2(1:L/2+1); 
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
%take natural log
P1 = log(P1);

plot(f,P1) 
title('Single-Sided Amplitude Spectrum of Signal')
xlabel('f (Hz)')
ylabel('Magnitude(log)')