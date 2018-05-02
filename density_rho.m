% Script for generating rho using CIC method

% Take floor for cell centers
i = floor(x);
j = floor(y);
k = floor(z);

% Now we need the distance from the actual particle to the cell center
dx = x - i ;
dy = y - j ;
dz = z - k ;

tx = 1 - dx;
ty = 1 - dy;
tz = 1 - dz;

rho = zeros(Ng,Ng,Ng);

for n = 1:Np
    
    ii = i(n)+1;
    jj = j(n)+1;
    kk = k(n)+1;
    
    % Special case 1, Particle at located at (Ng,Ng,Ng)
    if ii == Ng && jj == Ng && kk == Ng  
        
    rho(Ng,Ng,Ng) = ...
        rho(Ng,Ng,Ng) + mp*tx(n)*ty(n)*tz(n);
    
    rho(1,Ng,Ng) = ...
        rho(1,Ng,Ng) + mp*dx(n)*ty(n)*tz(n);
    
    rho(Ng,1,Ng) = ...
        rho(Ng,1,Ng) + mp*tx(n)*dy(n)*tz(n);
    
    rho(Ng,Ng,1) = ...
        rho(Ng,Ng,1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(1,1,Ng) = ...
        rho(1,1,Ng) + mp*dx(n)*dy(n)*tz(n);
    
    rho(1,Ng,1) = ...
        rho(1,Ng,1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(Ng,1,1) = ...
        rho(Ng,1,1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(1,1,1) = ...
        rho(1,1,1) + mp*dx(n)*dy(n)*dz(n);
    
    % Special case 2, Particle at located at (Ng,Ng,z)
    elseif ii == Ng && jj == Ng && kk < Ng
    
    rho(Ng,Ng,kk) = ...
        rho(Ng,Ng,kk) + mp*tx(n)*ty(n)*tz(n);
    
    rho(1,Ng,kk) = ...
        rho(1,Ng,kk) + mp*dx(n)*ty(n)*tz(n);
    
    rho(Ng,1,kk) = ...
        rho(Ng,1,kk) + mp*tx(n)*dy(n)*tz(n);
    
    rho(Ng,Ng,kk+1) = ...
        rho(Ng,Ng,kk+1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(1,1,kk) = ...
        rho(1,1,kk) + mp*dx(n)*dy(n)*tz(n);
    
    rho(1,Ng,kk+1) = ...
        rho(1,Ng,kk+1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(Ng,1,kk+1) = ...
        rho(Ng,1,kk+1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(1,1,kk+1) = ...
        rho(1,1,kk+1) + mp*dx(n)*dy(n)*dz(n);
    
    % Special case 3, Particle at located at (Ng,y,Ng)
    elseif ii == Ng && jj < Ng && kk == Ng 
    
    rho(Ng,jj,Ng) = ...
        rho(Ng,jj,Ng) + mp*tx(n)*ty(n)*tz(n);
    
    rho(1,jj,Ng) = ...
        rho(1,jj,Ng) + mp*dx(n)*ty(n)*tz(n);
    
    rho(Ng,jj+1,Ng) = ...
        rho(Ng,jj+1,Ng) + mp*tx(n)*dy(n)*tz(n);
    
    rho(Ng,jj,1) = ...
        rho(Ng,jj,1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(1,jj+1,Ng) = ...
        rho(1,jj+1,Ng) + mp*dx(n)*dy(n)*tz(n);
    
    rho(1,jj,1) = ...
        rho(1,jj,1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(Ng,jj+1,1) = ...
        rho(Ng,jj+1,1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(1,jj+1,1) = ...
        rho(1,jj+1,1) + mp*dx(n)*dy(n)*dz(n);
    
    % Special case 4, Particle at located at (x,Ng,Ng)
    elseif ii < Ng && jj == Ng && kk == Ng
    
    rho(ii,Ng,Ng) = ...
        rho(ii,Ng,Ng) + mp*tx(n)*ty(n)*tz(n);
    
    rho(ii+1,Ng,Ng) = ...
        rho(ii+1,Ng,Ng) + mp*dx(n)*ty(n)*tz(n);
    
    rho(ii,1,Ng) = ...
        rho(ii,1,Ng) + mp*tx(n)*dy(n)*tz(n);
    
    rho(ii,Ng,1) = ...
        rho(ii,Ng,1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(ii+1,1,Ng) = ...
        rho(ii+1,1,Ng) + mp*dx(n)*dy(n)*tz(n);
    
    rho(ii+1,Ng,1) = ...
        rho(ii+1,Ng,1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(ii,1,1) = ...
        rho(ii,1,1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(ii+1,1,1) = ...
        rho(ii+1,1,1) + mp*dx(n)*dy(n)*dz(n);
    
    % Special case 5, Particle at located at (x,y,Ng)
    elseif ii < Ng && jj< Ng && kk == Ng 
    
    rho(ii,jj,Ng) = ...
        rho(ii,jj,Ng) + mp*tx(n)*ty(n)*tz(n);
    
    rho(ii+1,jj,Ng) = ...
        rho(ii+1,jj,Ng) + mp*dx(n)*ty(n)*tz(n);
    
    rho(ii,jj+1,Ng) = ...
        rho(ii,jj+1,Ng) + mp*tx(n)*dy(n)*tz(n);
    
    rho(ii,jj,1) = ...
        rho(ii,jj,1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(ii+1,jj+1,Ng) = ...
        rho(ii+1,jj+1,Ng) + mp*dx(n)*dy(n)*tz(n);
    
    rho(ii+1,jj,1) = ...
        rho(ii+1,jj,1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(ii,jj+1,1) = ...
        rho(ii,jj+1,1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(ii+1,jj+1,1) = ...
        rho(ii+1,jj+1,1) + mp*dx(n)*dy(n)*dz(n);
        
    % Special case 6, Particle at located at (x,Ng,z)
    elseif ii < Ng && jj == Ng && kk < Ng
    
    rho(ii,Ng,kk) = ...
        rho(ii,Ng,kk) + mp*tx(n)*ty(n)*tz(n);
    
    rho(ii+1,Ng,kk) = ...
        rho(ii+1,Ng,kk) + mp*dx(n)*ty(n)*tz(n);
    
    rho(ii,1,kk) = ...
        rho(ii,1,kk) + mp*tx(n)*dy(n)*tz(n);
    
    rho(ii,Ng,kk+1) = ...
        rho(ii,Ng,kk+1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(ii+1,1,kk) = ...
        rho(ii+1,1,kk) + mp*dx(n)*dy(n)*tz(n);
    
    rho(ii+1,Ng,kk+1) = ...
        rho(ii+1,Ng,kk+1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(ii,1,kk+1) = ...
        rho(ii,1,kk+1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(ii+1,1,kk+1) = ...
        rho(ii+1,1,kk+1) + mp*dx(n)*dy(n)*dz(n);
    
    % Special case 6, Particle at located at (x,Ng,Ng)
    elseif ii ==  Ng && jj < Ng && kk < Ng
    
    rho(Ng,jj,kk) = ...
        rho(Ng,jj,kk) + mp*tx(n)*ty(n)*tz(n);
    
    rho(1,jj,kk) = ...
        rho(1,jj,kk) + mp*dx(n)*ty(n)*tz(n);
    
    rho(Ng,jj+1,kk) = ...
        rho(Ng,jj+1,kk) + mp*tx(n)*dy(n)*tz(n);
    
    rho(Ng,jj,kk+1) = ...
        rho(Ng,jj,kk+1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(1,jj+1,kk) = ...
        rho(1,jj+1,kk) + mp*dx(n)*dy(n)*tz(n);
    
    rho(1,jj,kk+1) = ...
        rho(1,jj,kk+1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(Ng,jj+1,kk+1) = ...
        rho(Ng,jj+1,kk+1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(1,jj+1,kk+1) = ...
        rho(1,jj+1,kk+1) + mp*dx(n)*dy(n)*dz(n);
    
    % All other cases finally
    
    else
    rho(ii,jj,kk) = ...
        rho(ii,jj,kk) + mp*tx(n)*ty(n)*tz(n);
    
    rho(ii+1,jj,kk) = ...
        rho(ii+1,jj,kk) + mp*dx(n)*ty(n)*tz(n);
    
    rho(ii,jj+1,kk) = ...
        rho(ii,jj+1,kk) + mp*tx(n)*dy(n)*tz(n);
    
    rho(ii,jj,kk+1) = ...
        rho(ii,jj,kk+1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(ii+1,jj+1,kk) = ...
        rho(ii+1,jj+1,kk) + mp*dx(n)*dy(n)*tz(n);
    
    rho(ii+1,jj,kk+1) = ...
        rho(ii+1,jj,kk+1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(ii,jj+1,kk+1) = ...
        rho(ii,jj+1,kk+1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(ii+1,jj+1,kk+1) = ...
        rho(ii+1,jj+1,kk+1) + mp*dx(n)*dy(n)*dz(n);
    
    end
end
