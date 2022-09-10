function [G_AR] = ARmodel(y, N, P)
G_AR = zeros(N-P, P);
% Iterate over rows
for i=1:N-P
    % Iterate over columns
    for j=1:P
        % Pass data into the matrix
        G_AR(i,j) = y(i-1+j);
    end
end
G_AR = fliplr(G_AR);
