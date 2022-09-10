function determine_window(signal)
%# Plot the different windowed response on a particular signal on the same 
%# graph to compare their performance and decide an optimum.
%determine length of sequence for transient note
N = size(signal,1);

%windows
w_hamm = signal.* hamming(N);
w_hann = signal.* hanning(N);
w_bartlett = signal.* bartlett(N);
w_kaiser = signal.* kaiser(N);

%use fvtool to visualise different forms of windows
fvtool(w_hamm, 1, w_hann, 2, w_bartlett, 3, w_kaiser, 4)