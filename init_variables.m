% Step 1 of the Particle Mesh 
% Initial variables and stuff

mp = 10;
Np = 1;
Ng = 64;

global OmegaL; global OmegaM;
OmegaL = 0.7; 
OmegaM = 0.3;

global x; global y; global z;
x = zeros(Np,1);
y = zeros(Np,1);
z = zeros(Np,1);

global vx; global vy; global vz;
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