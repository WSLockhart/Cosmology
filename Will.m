% evolution part of the program

% NOTE: We need to make all initial input variables GLOBAL

% Set number of particles and grid cells
Np = 1000;
Ng = 100;

global OmegaL; global OmegaM;
OmegaL = 1; 
OmegaM = 1;

% Declare position and velocity arrays for each particle
x = zeros(Np,1);
y = zeros(Np,1);
z = zeros(Np,1);

vx = zeros(Np,1);
vy = zeros(Np,1);
vz = zeros(Np,1);

phi = zeros(Ng,Ng,Ng);

%________________________________________________

a0 = 0.1;
stepsize = 0.01;
scalefactor = a0;

while scalefactor < 1   %start at a=a0, go until a=1 (today)

    %cycle through each particle to update velocities
    for n = 1:Ng %NOTE: need to mod for the edges!! 
       % calculate grad(phi) at its location (by averaging neighbors)
       GradPhiX = (phi(round(x(n))+1,round(y(n)),round(z(n))) - ...
           phi(round(x(n))-1,round(y(n)),round(z(n))))/2;
       GradPhiY = (phi(round(x(n)),round(y(n))+1,round(z(n))) - ...
           phi(round(x(n)),round(y(n))-1,round(z(n))))/2;
       GradPhiZ = (phi(round(x(n)),round(y(n)),round(z(n))+1) - ...
           phi(round(x(n)),round(y(n)),round(z(n))-1))/2; 

        % update velocities
        f = ff(scalefactor);
        vx(n) = vx(n) + f*GradPhiX*stepsize;
        vy(n) = vy(n) + f*GradPhiY*stepsize;
        vz(n) = vz(n) + f*GradPhiZ*stepsize;
    end
    
    %   update position using deltaX = a^-2 f(a)*p*deltaT 
    for n = 1:Ng %need to mod for the edges!! 
        x(n) = x(n) + (scalefactor^-2)*f*vx(n)*stepsize;
        y(n) = y(n) + (scalefactor^-2)*f*vy(n)*stepsize;
        z(n) = z(n) + (scalefactor^-2)*f*vz(n)*stepsize;
    end

    scalefactor = scalefactor + stepsize;
end


% the function 'ff' that relates the derivative of momentum to the gradient 
% of the potential. This is where our particular cosmology comes into play
function output = ff(a)
    global OmegaL; global OmegaM;
    output = ((1/a)*(OmegaM+OmegaL*a^3))^(1/2);
end













