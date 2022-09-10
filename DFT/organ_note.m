addpath('Audio_Week1');

%load organ.wav
[y, Fs] = audioread('organ.wav')
Nsamps = length(y);

%plot single-sided spectrum
signal_fft = fft(y);
L = length(signal_fft);

%only retaining one side of information
P2 = abs(signal_fft/L);
P1 = P2(1:L/2+1); 
P1(2:end-1) = 2*P1(2:end-1);
P1 = log(P1);
f = Fs*(0:(L/2))/L;

plot(f,P1) 
title('Single-Sided Amplitude Spectrum of Organ.wav')
xlim([0 22500])
xlabel('f (Hz)')
ylabel('Magnitude')


%find peaks and sort
[pks,locs] = findpeaks(P1, Fs,'MinPeakProminence', 0.00005, 'MinPeakDistance',0.001);
%find frequency components of the peaks
fValues = locs;

%clean the locs dataset by setting significant figures
locs = round(locs, 2, 'significant');
locs = sort(locs);


%iterate to find the frequencies
num_peaks = length(fValues);
%copy of frequency values
notes = locs; 
%initialise an empty array for identified notes
single_notes = [];

%iterate until all frequency components in notes identified
while length(notes) > 1
    %choose the lowest frequency component
    identified_freq = notes(1);
    %remove multiples of the lowest frequency component
    notes(~mod(notes, identified_freq)) = [];
    %add identified note into single_notes array
    single_notes =[single_notes, identified_freq];
end


%create duplicate of single_notes
fundamental_frequencies = single_notes;
num_notes = length(fundamental_frequencies)

%iterate for each element in fundamental_frequencies
for i = 1:num_notes 
    %for elements after the current element
    for j = i:num_notes
        lower_bound = single_notes(i) * 0.05;
        upper_bound = single_notes(i) * 0.95;
        %calculate the mod
        mod_value = mod(single_notes(j), single_notes(i));
        if (lower_bound < mod_value) && (mod_value < upper_bound)
            %delete from array if failing similarity criteria
        else
            fundamental_frequencies(j) = [];
            %decrement length by 1 after deletion
            num_notes = num_notes -1
        end
    end
end

%sort and remove duplicates
fundamental_frequencies = sort(unique(fundamental_frequencies));
fundamental_frequencies = [single_notes(1),fundamental_frequencies];
fundamental_frequencies = fundamental_frequencies * Fs

    
    