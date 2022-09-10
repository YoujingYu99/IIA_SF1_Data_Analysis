%create random noisy signal of length 10000
x_in=0.5*cos([1:10000]*pi/4)+sin([1:10000]*pi/100)+randn(1,10000);
%new array of same length as the signal, which is to be the output
y_out=0*x_in;

%set data length N and blocklength which is half of the data length
N=512;
overlap=256;

%partition x into data segments of length N and overlap length N/2
x=buffer(x_in,N,overlap);

%set number of samples and number of frames
[N_samps,N_frames]=size(x);

%create a matrix of N_frames columns.
%each column is a hanning window of length N so the matrix is N*N_frames
x_w=repmat(hanning(N),1,N_frames).*x;

%iterate over the each column until the second last
for frame_no=1:N_frames-2

    %fft on each column of hanning window
    X_w(:,frame_no)=fft(x_w(:,frame_no));
    
    %duplicate into a second matrix
    Y_w(:,frame_no)=X_w(:,frame_no);
    
    %attenuate to 0.1 for window frequencies from 2 to N/8 then copy into Y
    Y_w(2:N/8,frame_no)=0.1*X_w(2:N/8,frame_no);
    %attenuate to 0.2 for window frequencies from N/4+1 to N/2 then copy
    %into Y
    Y_w(N/4+1:N/2,frame_no)=0.2*X_w(N/4+1:N/2,frame_no);
    %calculate the conjugate of Y_w from 2 to N/2 rows
    %then fill N/2+2 to N with those values, in reverse order
    Y_w(N:-1:N/2+2,frame_no)=conj(Y_w(2:N/2,frame_no));
    %compute ifft of Y_w and copy into y_w
    y_w(:,frame_no)=ifft(Y_w(:,frame_no));
    
    %final result obtained by concatenating the frames and throwing away
    %the overlaps
    y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+y_w(:,frame_no)';
    
end    
    