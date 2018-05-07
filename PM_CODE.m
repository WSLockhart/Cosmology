% Particle Mesh Code to simulate large-scale structure formation in the universe


%{ 
    Next Steps:
        - fix coordinate inversion in contour plot
        - fix CIC method
        - convert to dimensionless variables and astronomical units
        * initial conditions!

        - what should a_init be? 

%}

clearvars;
close all;
% _______ 1) Initialize Variables _______ %

global particleMass; particleMass = 5;
%NOTE: Insert dimensionless variables & conversions to
%astronomical units 

global gridSize; gridSize = 64;
% particleNum = 32768;
% X = zeros(particleNum, 1);
% Y = zeros(particleNum, 1);
% Z = zeros(particleNum, 1); 
% Px = zeros(particleNum, 1);
% Py = zeros(particleNum, 1);
% Pz = zeros(particleNum, 1);

% _______ 2) Initial Conditions _______ %

a_init = 0.001;  %This corresponds roughly to z=1100, the time of the CMB (after matter-domination)
deltaA = 0.003;
global scalefactor; scalefactor = a_init;
global OmegaL; global OmegaM;
OmegaL = 0.7; 
OmegaM = 0.3;
global Phi; Phi = zeros(gridSize,gridSize,gridSize);
PhiK = zeros(gridSize,gridSize,gridSize);
Green = zeros(gridSize,gridSize,gridSize);  %Green function

%Single_Wave_Gen
Initial_Perturbations_3D %set initial X,Y,Z positions based on power spectrum


printmovie = true;
if printmovie
    vid = VideoWriter('sim.avi');
    open(vid);
end

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


% _______ 3) The Main Loop _______ %

numsteps = (1-a_init)/deltaA;
disp("running simulation...")
for framenum = 1:numsteps+1  %start at a=a0, go until a=1 (today)
    
    DENSITY_NGP
    POTENTIAL

    %update momenta
    for n = 1:particleNum 
        %NGP method:
        gradPhi = GRADS(mod(round(X(n)), gridSize), ...
                        mod(round(Y(n)), gridSize), ...
                        mod(round(Z(n)), gridSize));
        
%         % CIC method means we must average the gradient from all 8 neighbors 
%         
%         i = I(n); j = J(n); k = K(n);
%         % Now we need the distance from the actual particle to the cell center
%         dx = X(n)-I(n); dy = Y(n)-J(n); dz = Z(n)-K(n);
%         tx = 1-dx; ty = 1-dy; tz = 1-dz;
%         
%         gradPhi = GRADS(i,j,k)*tx*ty*tz + GRADS(mod(i+1,gridSize),j,k)*dx*ty*tz + ...
%                   GRADS(i,mod(j+1,gridSize),k)*tx*dy*tz + GRADS(i,j,mod(k+1,gridSize))*tx*ty*dz + ...
%                   GRADS(mod(i+1,gridSize),mod(j+1,gridSize),k)*dx*dy*tz + ...
%                   GRADS(mod(i+1,gridSize),j,mod(k+1,gridSize))*dx*ty*dz + ...
%                   GRADS(i,mod(j+1,gridSize),mod(k+1,gridSize))*tx*dy*dz + ...
%                   GRADS(mod(i+1,gridSize),mod(j+1,gridSize),mod(k+1,gridSize))*dx*dy*dz ;
        
        Px(n) = Px(n) - f(scalefactor)*gradPhi(1)*deltaA;  %update the velocities
        Py(n) = Py(n) - f(scalefactor)*gradPhi(2)*deltaA;
        Pz(n) = Pz(n) - f(scalefactor)*gradPhi(3)*deltaA;
    end
    
    % update positions 
    a_half = scalefactor + deltaA/2; % The "leapfrog" method uses the momenta 
    %from a(n-1/2) to update the positions at a(n)
    for n = 1:particleNum  
        X(n) = mod( X(n) + (1/a_half^2)*f(a_half)*Px(n)*deltaA, gridSize );
        Y(n) = mod( Y(n) + (1/a_half^2)*f(a_half)*Py(n)*deltaA, gridSize );
        Z(n) = mod( Z(n) + (1/a_half^2)*f(a_half)*Pz(n)*deltaA, gridSize );
    end
        
    % update plot
    scatter3(X,Y,Z,'.')
    axis([0 gridSize 0 gridSize 0 gridSize])
    ax=gca;
    mymovie(framenum) = getframe(ax);
    if printmovie
        writeVideo(vid,mymovie(framenum));
    end
    redshift = (1/scalefactor) - 1;
    disp("a = " + scalefactor + ", z = " + redshift)
%     disp(X)
    
    scalefactor = scalefactor + deltaA;
end

if printmovie
    implay(mymovie)
    close(vid);
end

fig2 = figure('Name','Potential');
    surf(Phi(:,:,gridSize/2))
    xlabel('x')
    ylabel('y')
disp("simulation complete.")
%---------------------------------------------------%


function output = f(a)
% the function 'ff' relates the derivative of momentum to the gradient 
% of the potential. This is where our particular cosmology comes into play
    global OmegaL; global OmegaM;
    output = ((1/a)*(OmegaM+OmegaL*a^3))^(-1/2);
end


function gradients = GRADS(i,j,k)
% returns the 3 components of the gradient of Phi at a given cell (i,j,k) 

  global gridSize; global Phi;
  % note the +1's because indicies start at 1
  gradients(1) = ( ...
      Phi( mod(i+1,gridSize)+1, j+1, k+1 ) ...
    - Phi( mod(i-1,gridSize)+1, j+1, k+1 ) )/2;
  gradients(2) = ( ...
      Phi( i+1, mod(j+1,gridSize)+1, k+1 ) ...
    - Phi( i+1, mod(j-1,gridSize)+1, k+1 ) )/2;
  gradients(3) = ( ...
      Phi( i+1, j+1, mod(k+1,gridSize)+1 ) ...
    - Phi( i+1, j+1, mod(k-1,gridSize)+1 ) )/2;  
  
end
