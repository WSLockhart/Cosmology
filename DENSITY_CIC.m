% Script for generating rho using CIC method

%clear the density array
Rho = zeros(gridSize,gridSize,gridSize);

% CIC Method: Take floor for cell centers
I = mod(floor(X), gridSize);
J = mod(floor(Y), gridSize);
K = mod(floor(Z), gridSize);

for n = 1:particleNum
    i = I(n); j = J(n); k = K(n);
    % Now we need the distance from the actual particle to the cell center
    dx = X(n)-I(n); dy = Y(n)-J(n); dz = Z(n)-K(n);
    tx = 1-dx; ty = 1-dy; tz = 1-dz;

    %NOTE: The COORDINATES x ranges from 0 to N, so round(x) returns 0 to N-1
    %But the ARRAY Rho[] has indicies from 1 to N, so Rho[x+1] corresponds to location x.
    Rho(i+1, j+1, k+1) = ...
    Rho(i+1, j+1, k+1) + particleMass * tx * ty * tz; 
    
    Rho(mod(i+1,gridSize)+1, j+1, k+1) = ...
    Rho(mod(i+1,gridSize)+1, j+1, k+1) + particleMass * dx * ty * tz; 
    
    Rho(i+1, mod(j+1,gridSize)+1, k+1) = ...
    Rho(i+1, mod(j+1,gridSize)+1, k+1) + particleMass * tx * dy * tz; 
    
    Rho(i+1, j+1, mod(k+1,gridSize)+1) = ...
    Rho(i+1, j+1, mod(k+1,gridSize)+1) + particleMass * tx * ty * dz; 
    
    Rho(mod(i+1,gridSize)+1, mod(j+1,gridSize)+1, k+1) = ...
    Rho(mod(i+1,gridSize)+1, mod(j+1,gridSize)+1, k+1) + particleMass * dx * dy * tz; 
    
    Rho(mod(i+1,gridSize)+1, j+1, mod(k+1,gridSize)+1) = ...
    Rho(mod(i+1,gridSize)+1, j+1, mod(k+1,gridSize)+1) + particleMass * dx * ty * dz; 
    
    Rho(i+1, mod(j+1,gridSize)+1, mod(k+1,gridSize)+1) = ...
    Rho(i+1, mod(j+1,gridSize)+1, mod(k+1,gridSize)+1) + particleMass * tx * dy * dz; 
   
    Rho(mod(i+1,gridSize)+1, mod(j+1,gridSize)+1, mod(k+1,gridSize)+1) = ...
    Rho(mod(i+1,gridSize)+1, mod(j+1,gridSize)+1, mod(k+1,gridSize)+1) + particleMass * dx * dy * dz;

end









