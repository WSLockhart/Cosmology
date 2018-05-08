clearvars;

global gridSize; gridSize = 64;
obj_side = 32;
particleNum = obj_side ^ 3;
particleMass = 5; %MAKE RIGOROUS

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
Pz = zeros(particleNum, 1);