clearvars;

global gridSize; gridSize = 64;
obj_side = 32;
particleNum = obj_side ^ 3;

a_init = 0.001;  %This corresponds roughly to z=1100, the time of the CMB (after matter-domination)
global scalefactor; scalefactor = a_init;

global OmegaL; global OmegaM;
OmegaL = 0.7; 
OmegaM = 0.3;
global Phi; Phi = zeros(gridSize,gridSize,gridSize);
PhiK = zeros(gridSize,gridSize,gridSize);
Green = zeros(gridSize,gridSize,gridSize);  %Green function

X = zeros(particleNum, 1);
Y = zeros(particleNum, 1);
Z = zeros(particleNum, 1); 
Px = zeros(particleNum, 1);
Py = zeros(particleNum, 1);
Pz = zeros(particleNum, 1

% Units and stuff
% Constants:
Rho_crit = 2.735e41;    % Critical mass density in kg / Mpc^3
Hub_time = 4.55e17;     % Hubble time in seconds (time scale)

% Parameters:
cell_side_length = 1;   % Cell side length of 1 Mpc, r0

% Derived 
Rho_0 = Rho_crit*OmegaM;    % Characteristic mass density kg / Mpc^3
particleMass = 8*Rho_0;     % Mass of particles in kg

V_0 = cell_side_length / Hub_time;  % characteristic velocity
Phi_0 = V_0^2;
