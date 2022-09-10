%set length to be investigated
length_array = linspace(256, 256*4, 3);
length_of_array = length(length_array);
normalised_freq = 0.5;
N = 32;
n = 1:N;


%oscillating function
w_cos = cos(normalised_freq*pi*n);
%windows
w_hamm = w_cos' .* hamming(N);
w_hann = w_cos' .* hanning(N);
w_bartlett = w_cos' .* bartlett(N);
w_kaiser = w_cos' .* kaiser(N);


fvtool(w_cos,1,w_hamm,2,w_hann,3,w_bartlett,4,w_kaiser,5);

%choose one window function to analyse how it alters with length N
length_array = linspace(256, 256*2, 2+1);
length_of_array = length(length_array);

%find peaks
cos_peak = zeros(3);
hamm_peak = zeros(3);
hann_peak = zeros(3);
bartlett_peak = zeros(3);
kaiser_peak = zeros(3);
hamm_1_peak = zeros(3);
hamm_2_peak = zeros(3);

for i = 1:length_of_array
    N_chosen = length_array(i);
    n = 1:N_chosen;
    
    %oscillating function
    w_cos = cos(normalised_freq*pi*n);
    %windows
    w_hamm = w_cos' .* hamming(N_chosen);
    w_hann = w_cos' .* hanning(N_chosen);
    w_bartlett = w_cos' .* bartlett(N_chosen);
    w_kaiser = w_cos' .* kaiser(N_chosen);
    fv_plot = fvtool(w_cos, 1 ,w_hamm,2,w_hann,3,w_bartlett,4,w_kaiser,5);
    legend(fv_plot, 'original', 'Hamming', 'Hanning', 'Bartlett', 'Kaiser');
    
    %find peaks
    [pks_hamm,locs_hamm] = findpeaks(w_hamm, normalised_freq*pi*n, 'MinPeakHeight', 0.9999);
    [pks_hann,locs_hann] = findpeaks(w_hann, normalised_freq*pi*n, 'MinPeakHeight', 0.9999);
    [pks_bartlett,locs_bartlett] = findpeaks(w_bartlett, normalised_freq*pi*n, 'MinPeakHeight', 0.99);
    [pks_kaiser,locs_kaiser] = findpeaks(w_kaiser, normalised_freq*pi*n, 'MinPeakHeight', 0.99999);
    hamm_peak(i) = locs_hamm / (N_chosen * normalised_freq);
    hann_peak(i) = locs_hann/ (N_chosen * normalised_freq);
    bartlett_peak(i) = locs_bartlett/ (N_chosen * normalised_freq);
    kaiser_peak(i) = locs_kaiser/ (N_chosen * normalised_freq);
  

    %add noise
    n1 = 0.5*randn(N_chosen,1)';
    n2 = 1*randn(N_chosen,1)';
    w_cos_N_1 = cos(normalised_freq*pi*n) + n1;
    w_hamm_1 = w_cos_N_1'.* hamming(N_chosen);
    w_cos_N_2 = cos(normalised_freq*pi*n) + n2;
    w_hamm_2 = w_cos_N_2'.* hamming(N_chosen);
    w_cos_control = cos(normalised_freq*pi*n);
    w_hamm_control = w_cos_control'.* hamming(N_chosen);
    fv_noise = fvtool(w_hamm_control, 1, w_hamm_1, 2, w_hamm_2, 3);
    legend(fv_noise, 'Original', 'Standard Deviation 0.5', 'Standard Deviation 1');
    
    %find peaks
    [pks_hamm,locs_hamm_1] = findpeaks(w_hamm_1, normalised_freq*pi*n, 'MinPeakHeight', 1.6);
    [pks_hamm,locs_hamm_2] = findpeaks(w_hamm_1, normalised_freq*pi*n, 'MinPeakHeight', 1.6);
    hamm_1_peak(i) = locs_hamm_1 / (N_chosen * normalised_freq);
    hamm_2_peak(i) = locs_hamm_2/ (N_chosen * normalised_freq);
  
end


