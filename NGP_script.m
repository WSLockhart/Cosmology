% Particle Mesh code for Astro class by MEEEE Ed Schenk  Mar. 30 2018

% Set number of particles and grid cells
Np = 8;
Ng = 128;

% Set some variables
mp = 1;
Omega = 1;
OmegaM = 1;
a = 1;
H0 = 1;
G = 1;
rho0 = 1;
%rho0 = (3*H0^2/(8*pi*G))*OmegaM;

% Declare position and velocity arrays for each particle
x = zeros(Np,1);
y = zeros(Np,1);
z = zeros(Np,1);

vx = zeros(Np,1);
vy = zeros(Np,1);
vz = zeros(Np,1);

% Declaring density and potential arrays for each grid point
rho = zeros(Ng,Ng,Ng);
phi = zeros(Ng,Ng,Ng);
phik= zeros(Ng,Ng,Ng);
% Declaring Green's function
Gr = zeros(Ng,Ng,Ng);

% Let's use Nearest Grid Point (NGP) method
% setting random particle positons

x = Ng*rand(Np,1);
y = Ng*rand(Np,1);
z = Ng*rand(Np,1);

% Plotting initial points
%scatter3(x,y,z)
axis([0 Ng 0 Ng 0 Ng])
hold on

% Round to get nearest grid points
i = mod(round(x), Ng);
j = mod(round(y), Ng);
k = mod(round(z), Ng);

% Placing each point into the density array
for n = 1:Np
    rho(i(n)+1,j(n)+1,k(n)+1)= rho(i(n)+1,j(n)+1,k(n)+1) + 1;
end

% Average density
rhoav = mean(rho(:));

deltax = rho - 1;

% Plotting nearest grid points
scatter3(i,j,k);
hold off

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
contourslice(phi,i,j,k)
view(3)
grid on

