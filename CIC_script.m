% Script for CIC (cloud in cell) method By Ed Apr. 12, 2017

mp = 1;
Np = 24;
Ng = 128;

x = zeros(Np,1);
y = zeros(Np,1);
z = zeros(Np,1);

rho = zeros(Ng,Ng,Ng);

x = Ng*rand(Np,1);
y = Ng*rand(Np,1);
z = Ng*rand(Np,1);

% New stuff

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
    
    rho(mod(i(n)+1,Ng)+1,j(n)+1,k(n)+1) = ...
        rho(mod(i(n)+1,Ng)+1,j(n)+1,k(n)+1) + mp*dx(n)*ty(n)*tz(n);
    
    rho(i(n)+1,mod(j(n)+1,Ng)+1,k(n)+1) = ...
        rho(i(n)+1,mod(j(n)+1,Ng)+1,k(n)+1) + mp*tx(n)*dy(n)*tz(n);
    
    rho(i(n)+1,j(n)+1,mod(k(n)+1,Ng)+1) = ...
        rho(i(n)+1,j(n)+1,mod(k(n)+1,Ng)+1) + mp*tx(n)*ty(n)*dz(n);
    
    rho(mod(i(n)+1,Ng)+1,mod(j(n)+1,Ng)+1,k(n)+1) = ...
        rho(mod(i(n)+1,Ng)+1,mod(j(n)+1,Ng)+1,k(n)+1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(mod(i(n)+1,Ng)+1,j(n)+1,mod(k(n)+1,Ng)+1) = ...
        rho(mod(i(n)+1,Ng)+1,j(n)+1,mod(k(n)+1,Ng)+1) + mp*dx(n)*ty(n)*dz(n);
    
    rho(mod(i(n)+1,Ng)+1,mod(j(n)+1,Ng)+1,k(n)+1) = ...
        rho(mod(i(n)+1,Ng)+1,mod(j(n)+1,Ng)+1,k(n)+1) + mp*dx(n)*dy(n)*tz(n);
    
    rho(mod(i(n)+1,Ng)+1,mod(j(n)+1,Ng)+1,mod(k(n)+1,Ng)+1) = ...
        rho(mod(i(n)+1,Ng)+1,mod(j(n)+1,Ng)+1,mod(k(n)+1,Ng)+1) + mp*dx(n)*dy(n)*dz(n);
    
end


% CIC complete
% Previous code pasted below


phik= zeros(Ng,Ng,Ng);
Gr = zeros(Ng,Ng,Ng);

% Average density
rhoav = mean(rho(:));

deltax = rho - 1;

% Start solving the discretized Poisson Equation
% Need to solve in k-space
% Get the overdensity in k-space
deltak = fftn(deltax);

% Determine green's function array
for l = 1:Ng
    for m = 1:Ng
        for n = 1:Ng
            if (l==1)&&(m==1)&&(n==1)
                Gr(1,1,1)=0;
            else
                Gr(l,m,n)= -(1/(sin(pi*(l-1)/Ng)^2+sin(pi*(m-1)/Ng)^2+sin(pi*(n-1)/Ng)^2));
            end
        end
    end
end

% Finding phi in fourier space
for l = 1:Ng
    for m = 1:Ng
        for n = 1:Ng
            phik(m,l,n)= Gr(l,m,n)*deltak(l,m,n);
        end
    end
end

% Now we FFT phik back to configuration space to get phi
phi = circshift(circshift(circshift(...  Shifting is to align the potential field to the mass positions
    real( ...   potential must be real
    ifftn(phik) ...
    ),-1,1),-1,2),-1,3);

%plotting bullshit
hold off
contourslice(phi,x,y,z)
hold on
scatter3(x,y,z,'filled')
axis([0 Ng 0 Ng 0 Ng])
grid on
