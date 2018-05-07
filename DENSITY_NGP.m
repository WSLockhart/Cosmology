% Script for generating rho using NGP method
global particleMass;

%clear the density array
Rho = zeros(gridSize,gridSize,gridSize);

Xav = mod(round(X), gridSize);
Yav = mod(round(Y), gridSize);
Zav = mod(round(Z), gridSize);       

for n = 1:particleNum
    %NOTE: The COORDINATES x ranges from 0 to N, so round(x) returns 0 to N-1
    %But the ARRAY Rho[] has indicies from 1 to N, so Rho[x+1] corresponds
    %to location x.
    Rho(Xav(n)+1,Yav(n)+1,Zav(n)+1) = ...
    Rho(Xav(n)+1,Yav(n)+1,Zav(n)+1) + particleMass;  %simple NGP method
end
