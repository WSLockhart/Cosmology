% Script for generating rho using CIC method

%clear the density array
Rho = zeros(gridSize,gridSize,gridSize);

% Take floor for cell centers
I = mod(floor(X), gridSize);
J = mod(floor(Y), gridSize);
K = mod(floor(Z), gridSize);

% Now we need the distance from the actual particle to the cell center
dx = X - I;
dy = Y - J;
dz = Z - K;

tx = 1 - dx;
ty = 1 - dy;
tz = 1 - dz;

for n = 1:particleNum
    %NOTE: The COORDINATES x ranges from 0 to N, so round(x) returns 0 to N-1
    %But the ARRAY Rho[] has indicies from 1 to N, so Rho[x+1] corresponds to location x.
    Rho(I(n)+1, J(n)+1, K(n)+1) = ...
    Rho(I(n)+1, J(n)+1, K(n)+1) + particleMass * tx * ty * tz; 
    
    Rho(mod(I(n)+1,gridSize)+1, J(n)+1, K(n)+1) = ...
    Rho(mod(I(n)+1,gridSize)+1, J(n)+1, K(n)+1) + particleMass * dx * ty * tz; 
    
    Rho(I(n)+1, mod(J(n)+1,gridSize)+1, K(n)+1) = ...
    Rho(I(n)+1, mod(J(n)+1,gridSize)+1, K(n)+1) + particleMass * tx * dy * tz; 
    
    Rho(I(n)+1, J(n)+1, mod(K(n)+1,gridSize)+1) = ...
    Rho(I(n)+1, J(n)+1, mod(K(n)+1,gridSize)+1) + particleMass * tx * ty * dz; 
    
    Rho(mod(I(n)+1,gridSize)+1, mod(J(n)+1,gridSize)+1, K(n)+1) = ...
    Rho(mod(I(n)+1,gridSize)+1, mod(J(n)+1,gridSize)+1, K(n)+1) + particleMass * dx * dy * tz; 
    
    Rho(mod(I(n)+1,gridSize)+1, J(n)+1, mod(K(n)+1,gridSize)+1) = ...
    Rho(mod(I(n)+1,gridSize)+1, J(n)+1, mod(K(n)+1,gridSize)+1) + particleMass * dx * ty * dz; 
    
    Rho(I(n)+1, mod(J(n)+1,gridSize)+1, mod(K(n)+1,gridSize)+1) = ...
    Rho(I(n)+1, mod(J(n)+1,gridSize)+1, mod(K(n)+1,gridSize)+1) + particleMass * tx * dy * dz; 
   
    Rho(mod(I(n)+1,gridSize)+1, mod(J(n)+1,gridSize)+1, mod(K(n)+1,gridSize)+1) = ...
    Rho(mod(I(n)+1,gridSize)+1, mod(J(n)+1,gridSize)+1, mod(K(n)+1,gridSize)+1) + particleMass * dx * dy * dz;

end









