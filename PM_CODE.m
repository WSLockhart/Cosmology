% NEW SCRIPT
% We are starting our coordinates at zero for god's sake

%{ 
    Next Steps:
        - fix coordinate inversion in contour plot
        - switch to CIC method over NGP
        - re-introduce the scale factor into the equations
        - convert to dimensionless variables and astronomical units
        * initial conditions!

    (minor):
        - close figures between runs
%}

clearvars
% _______ 1) Initialize Variables _______ %
particleMass = 10;
particleNum = 1;
gridSize = 64; %NOTE: Insert dimensionless variables & conversions to
%astronomical units 

printmovie = false;

global OmegaL; global OmegaM;
OmegaL = 0.7; 
OmegaM = 0.3;

% _______ 2) Initial Conditions _______ %

% (random for now - this will be replaced by an initial power spectrum)
X = gridSize*rand(particleNum,1);
Y = gridSize*rand(particleNum,1);   %pNx1 array of random numbers between 0 and gridSize
Z = zeros(particleNum,1); Z = Z + gridSize/2;
% X = [20.3];
% Y = [gridSize/2];
% Z = [gridSize/2];

Vx = zeros(particleNum,1);
Vy = zeros(particleNum,1);
Vz = zeros(particleNum,1);

Phi = zeros(gridSize,gridSize,gridSize);
PhiK = zeros(gridSize,gridSize,gridSize);
Green = zeros(gridSize,gridSize,gridSize);  %Green function

DENSITY_NGP
POTENTIAL

%Pl0TS:
fig1 = figure('Name','INITIAL DATA');
    subplot(2,1,1) % add first plot in 2 x 1 grid
    contourf(Phi(:,:,gridSize/2))
    title('Potential Contours');  
    subplot(2,1,2) % add second plot in 2 x 1 grid
    plot1 = scatter(X,Y);
    axis([0 gridSize 0 gridSize])
    title('Particles');
    set(gcf, 'Position', [1200, 1000, 400, 800])

fig2 = figure('Name','Potential');
    surf(Phi(:,:,gridSize/2))
    
% _______ 3) The Main Loop _______ %
a0 = 1; 
stepsize = 0.01;
scalefactor = a0;
framenum = 1;

disp("running simulation...")
while scalefactor < 1   %start at a=a0, go until a=1 (today)
    
    DENSITY_NGP
    POTENTIAL

    % for the NGP method we just need to round - update to CIC later
    I = mod(round(X), gridSize);
    J = mod(round(Y), gridSize);
    K = mod(round(Z), gridSize);
    
    %cycle through each particle and update velocities
    for n = 1:particleNum
        gradPhiX = ( ...
            Phi( mod(I(n)+1,gridSize)+1, J(n)+1, K(n)+1 ) ...
          - Phi( mod(I(n)-1,gridSize)+1, J(n)+1, K(n)+1 ) )/2;
      % note the +1's because indicies start at 1
         gradPhiY = ( ...
            Phi( I(n)+1, mod(J(n)+1,gridSize)+1, K(n)+1 ) ...
          - Phi( I(n)+1, mod(J(n)-1,gridSize)+1, K(n)+1 ) )/2;
        gradPhiZ=0;
        
        f = ff(scalefactor);
        Vx(n) = Vx(n) - f*gradPhiX*stepsize;
        Vy(n) = Vy(n) - f*gradPhiY*stepsize;
        Vz(n) = Vz(n) - f*gradPhiZ*stepsize;
    end
    
    % update positions 
    for n = 1:particleNum  
%       X(n) = X(n) + (scalefactor^-2)*f*vx(n)*stepsize;
        X(n) = mod( X(n) + Vx(n)*stepsize,gridSize );
        Y(n) = mod( Y(n) + Vy(n)*stepsize,gridSize );
        Z(n) = mod( Z(n) + Vz(n)*stepsize,gridSize );
    end
    
    % update positions plot
    scatter(X,Y)
    axis([0 gridSize 0 gridSize])
    ax=gca;
    mymovie(framenum) = getframe(ax);
   
    
    scalefactor = scalefactor + stepsize;
    framenum = framenum + 1;
    disp("a = " + scalefactor)
end

if printmovie
    implay(mymovie)
end

disp("simulation complete.")


% the function 'ff' relates the derivative of momentum to the gradient 
% of the potential. This is where our particular cosmology comes into play
function output = ff(a)
    global OmegaL; global OmegaM;
%   output = ((1/a)*(OmegaM+OmegaL*a^3))^(1/2);
    output = 0.2;  
end