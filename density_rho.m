% Script for generating rho using CIC method

% Take floor for cell centers
i = floor(x);
j = floor(y);
k = floor(z);

% Now we need the distance from the actual particle to the cell center
dx = x - i;
dy = y - j;
dz = z - k;

tx = 1 - dx;
ty = 1 - dy;
tz = 1 - dz;

for n = 1:Np
    rho(i(n)+1,j(n)+1,k(n)+1) = ...
        rho(i(n)+1,j(n)+1,k(n)+1) + mp*tx(n)*ty(n)*tz(n);
    
    rho(mod(i(n)+1,Ng+1)+1,j(n)+1,k(n)+1) = ...
        rho(mod(i(n)+1,Ng+1)+1,j(n)+1,k(n)+1) + mp*dx(n)*ty(n)*tz(n);
    
    rho(i(n)+1,mod(j(n)+1,Ng+1)+1,k(n)+1) = ...
        rho(i(n)+1,mod(j(n)+1,Ng+1)+1,k(n)+1) + mp*tx(n)*dy(n)*tz(n);
    
    rho(i(n)+1,j(n)+1,mod(k(n)+1,Ng+1)+1) = ...
        rho(i(n)+1,j(n)+1,mod(k(n)+1,Ng+1)+1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(mod(i(n)+1,Ng+1)+1,mod(j(n)+1,Ng+1)+1,k(n)+1) = ...
        rho(mod(i(n)+1,Ng+1)+1,mod(j(n)+1,Ng+1)+1,k(n)+1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(mod(i(n)+1,Ng+1)+1,j(n)+1,mod(k(n)+1,Ng+1)+1) = ...
        rho(mod(i(n)+1,Ng+1)+1,j(n)+1,mod(k(n)+1,Ng+1)+1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(mod(i(n)+1,Ng+1)+1,mod(j(n)+1,Ng+1)+1,k(n)+1) = ...
        rho(mod(i(n)+1,Ng+1)+1,mod(j(n)+1,Ng+1)+1,k(n)+1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(mod(i(n)+1,Ng+1)+1,mod(j(n)+1,Ng+1)+1,mod(k(n)+1,Ng+1)+1) = ...
        rho(mod(i(n)+1,Ng+1)+1,mod(j(n)+1,Ng+1)+1,mod(k(n)+1,Ng+1)+1) + mp*dx(n)*dy(n)*dz(n);
    
end
