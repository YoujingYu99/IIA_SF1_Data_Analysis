function X = DFT(signal,N_chosen)
N=N_chosen;
X=zeros(1,N);
for p = 1:N
    X_p = 0;
    for n = 1:N
        x_indiv=exp(-1j*2*pi*(p-1)*(n-1)/N)*signal(n);
        X_p = X_p + x_indiv;
    end
X(p)=X_p;
end
end
    

