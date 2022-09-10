function X = DFT_vectorized(signal,N_chosen)
%# Calculate DFT of input signal using meshgrid mapping with the datalength
%# chosen.
%# INPUT signal to perform DFT on
%# OUTPUT DFT of signal
%offsetting ranges of p and n
p = 0:N_chosen - 1;
n = 0:N_chosen - 1;
%creating meshgrid of p and n
[p2d, n2d] = meshgrid(p, n);
DFT_mat = exp(-2j*pi*p2d.*n2d/N_chosen);
%applying the DFT matrix on the signal
X = DFT_mat*signal(1:N_chosen);
end
    