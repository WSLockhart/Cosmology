% The Master Script 4/16/2017 V0.01

% Run step1.m
% Run random_masses.m
% Run density_rho.m
% Run potential_phi.m
% Run Will.m

% 1)
init_variables

% 2)
%random_masses
random_masses2D
% scatter(x,y)
% axis([0 Ng 0 Ng])

% Loop these sections

a0 = 0.1;
stepsize = 0.05;
scalefactor = a0;
framenum = 1;

while scalefactor < 1   %start at a=a0, go until a=1 (today)
   
    % 3)
    density_rho

    % 4)potential_phi
    potential_phi2D    

    %cycle through each particle to update velocities
    for n = 1:Np
        
    % 1. calculate grad(phi) at its location (by averaging neighbors)
       ii = i(n);
       jj = j(n);
       kk=1;  % JERRY-RIG FOR 2D
       if ii == 1
           GradPhiX = ( phi(ii+1,jj,kk) - phi(Ng,jj,kk) )/2;
       elseif ii == Ng
           GradPhiX = ( phi(1,jj,kk) - phi(ii-1,jj,kk) )/2;
       else
           GradPhiX = ( phi(ii+1,jj,kk) - phi(ii-1,jj,kk) )/2;
       end
   % NOTE: SUPER GHETTO - USE MOD INSTEAD
       if jj == 1
           GradPhiY = ( phi(ii,jj+1,kk) - phi(ii,Ng,kk) )/2;
       elseif jj == Ng
           GradPhiY = ( phi(ii,1,kk) - phi(ii,jj-1,kk) )/2;
       else
           GradPhiY = ( phi(ii,jj+1,kk) - phi(ii,jj-1,kk) )/2;
       end
   % NOTE: SUPER GHETTO - USE MOD INSTEAD
     
       
    % 3. update velocities
       f = ff(scalefactor);
       vx(n) = vx(n) - f*GradPhiX*stepsize;
       vy(n) = vy(n) - f*GradPhiY*stepsize;
    end
    
% 3. update position using deltaX = a^-2 f(a)*p*deltaT 
    for n = 1:Np  
        x(n) = x(n) + (scalefactor^-2)*f*vx(n)*stepsize;
        if x(n) < 1
            x(n) = Ng-1+x(n);
        elseif x(n) > Ng
            x(n) = x(n)-Ng+1;
        end
      % NOTE: SUPER GHETTO - USE MOD INSTEAD
        y(n) = y(n) + (scalefactor^-2)*f*vy(n)*stepsize;
        if y(n) < 1
            y(n) = Ng-1+y(n);
        elseif y(n) > Ng
            y(n) = y(n)-Ng+1;
        end
        
    end

    
    %     %plotting bullshit
%     hold off
%     contourf(phi(:,:,1))
%     hold on
%     %scatter3(x,y,z)
%     axis([0 Ng 0 Ng 0 Ng])
%     grid on
%  
    scatter(x,y)
    axis([1 Ng 1 Ng])
    ax=gca;
    mymovie(framenum) = getframe(ax);


    scalefactor = scalefactor + stepsize;
    framenum = framenum + 1;
end

figure
movie(mymovie,5)
    
% the function 'ff' that relates the derivative of momentum to the gradient 
% of the potential. This is where our particular cosmology comes into play
function output = ff(a)
    global OmegaL; global OmegaM;
    output = ((1/a)*(OmegaM+OmegaL*a^3))^(1/2);
end
