% Script for generating rho using NGP method

%clear the density array
Rho = zeros(gridSize,gridSize,gridSize);

% Round to get nearest grid points
I = mod(round(X), gridSize);
J = mod(round(Y), gridSize);
K = mod(round(Z), gridSize);

for n = 1:particleNum
    %NOTE: The COORDINATES x ranges from 0 to N, so round(x) returns 0 to N-1
    %But the ARRAY Rho[] has indicies from 1 to N, so Rho[x+1] corresponds
    %to location x.
    Rho(I(n)+1,J(n)+1,K(n)+1) = Rho(I(n)+1,J(n)+1,K(n)+1) + particleMass;  %simple NGP method
end
