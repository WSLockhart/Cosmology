% The Master Script 4/16/2017 V0.01

% Run step1.m
% Run random_masses.m
% Run density_rho.m
% Run potential_phi.m
% Run Will.m
global x,global y,global z;
global vx; global vy; global vz;

clearvars

% 1)
init_variables

% 2)
random_masses2D


density_rho
potential_phi2D   

figure('Name','INITIAL')
    subplot(2,1,1)       % add first plot in 2 x 1 grid
    contourf(phi(:,:,1))
    title('Potential Contours');  
    subplot(2,1,2)       % add second plot in 2 x 1 grid
    plot1 = scatter(x,y);
    axis([0 Ng 0 Ng])
    title('Particles');


% 3) The Main Loop
a0 = 0.1;
stepsize = 0.05;
scalefactor = a0;
framenum = 1;

while scalefactor < 3   %start at a=a0, go until a=1 (today)
%     disp(x)

    % calculate the density of each cell using the CIC method
    density_rho

    % calculate the potential everywhere
    potential_phi2D    

    %cycle through each particle to update velocities
    for n = 1:Np
        
    % 1. calculate grad(phi) at its location (by averaging neighbors)
       ii = i(n);
       jj = Ng/2;
       kk = Ng/2;  % JERRY-RIG FOR 2D
       if ii == 1
           GradPhiX = ( phi(ii+1,jj,kk) - phi(Ng,jj,kk) )/2;
       elseif ii == Ng
           GradPhiX = ( phi(1,jj,kk) - phi(ii-1,jj,kk) )/2;
       else
           GradPhiX = ( phi(ii+1,jj,kk) - phi(ii-1,jj,kk) )/2;
       end
   % NOTE: SUPER GHETTO - USE MOD INSTEAD
%        if jj == 1
%            GradPhiY = ( phi(ii,jj+1,kk) - phi(ii,Ng,kk) )/2;
%        elseif jj == Ng
%            GradPhiY = ( phi(ii,1,kk) - phi(ii,jj-1,kk) )/2;
%        else
%            GradPhiY = ( phi(ii,jj+1,kk) - phi(ii,jj-1,kk) )/2;
%        end
  
     
%    figure('Name','Grad Phi');  
%    quiver(x,y,GradPhiX,GradPhiY)
   
    % 3. update velocities
%        f = ff(scalefactor);
       f=1;
       vx(n) = vx(n) - f*GradPhiX*stepsize;
%        vy(n) = vy(n) - f*GradPhiY*stepsize;
    end

% 3. update position using deltaX = a^-2 f(a)*p*deltaT 
    for n = 1:Np  
%         x(n) = x(n) + (scalefactor^-2)*f*vx(n)*stepsize*100;
        x(n) = x(n) + vx(n)*stepsize*100;
        if x(n) < 1
            x(n) = Ng-1+x(n);
        elseif x(n) > Ng
            x(n) = x(n)-Ng+1;
        end
%       % NOTE: SUPER GHETTO - USE MOD INSTEAD
%         y(n) = y(n) + (scalefactor^-2)*f*vy(n)*stepsize*10;
%         if y(n) < 1
%             y(n) = Ng-1+y(n);
%         elseif y(n) > Ng
%             y(n) = y(n)-Ng+1;
%         end   
    end

    scatter(x,y)
    axis([0 Ng 0 Ng])
    ax=gca;
    mymovie(framenum) = getframe(ax);

    scalefactor = scalefactor + stepsize;
    framenum = framenum + 1;
end


implay(mymovie)

% 
% figure('Name','FINAL')
%     subplot(2,1,1)       % add first plot in 2 x 1 grid
%     plot1 = scatter(x,y);
%     axis([0 Ng 0 Ng])
%     title('Particles');
%     subplot(2,1,2)       % add second plot in 2 x 1 grid
%     contourf(phi(:,:,1))
%     title('Potential Contours');
% 

% the function 'ff' that relates the derivative of momentum to the gradient 
% of the potential. This is where our particular cosmology comes into play
function output = ff(a)
    global OmegaL; global OmegaM;
    output = ((1/a)*(OmegaM+OmegaL*a^3))^(1/2);
end


