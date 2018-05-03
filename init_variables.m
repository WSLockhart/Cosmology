% Step 1 of the Particle Mesh 
% Initial variables and stuff

global OmegaL; global OmegaM;
OmegaL = 0.7; 
OmegaM = 0.3;

% Units:
% Length in Mpc
% mass in solar masses
% time in seconds
% velocity in Mpc/s

% mass of each galaxy in solar masses
mp = 1e12;

% Gravitational constant G in units of Mpc^3/s^2/Solarmass
grav_const = 4.515e-48;

% Length of the box in mega parsecs
Length = 100;

% Unit length (length of a grid cell in Mpc)
r0 = Length/Ng;

% Unit of time (Hubble time in seconds)
t0 = 4.55e17;

% Unit of velocity (in units Mpc/s)
v0 = r0/t0;

% Unit of density (units of solar masses per cubic Mpc)
rho0 = (3*OmegaM)/(8*pi*t0^2*grav_const);

% Unit of potential phi (units of Mpc^2/s^2)
phi0 = v0^2;

% Number of grid points and number of particles
Np = 2;
Ng = 64;

x = zeros(Np,1);
y = zeros(Np,1);
z = zeros(Np,1);

i = zeros(Np,1);
j = zeros(Np,1);
k = zeros(Np,1);

vx = zeros(Np,1);
vy = zeros(Np,1);
vz = zeros(Np,1);

rho = zeros(Ng,Ng,Ng);
phi = zeros(Ng,Ng,Ng);
Gr  = zeros(Ng,Ng,Ng);

% Determine green's function array Gr, for calculating phi
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

