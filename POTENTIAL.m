% Calculate the gravitational potential everywhere

% RhoAv = mean(Rho(:));
% DeltaX = (Rho - RhoAv)/RhoAv; % over-density in position space
DeltaX = Rho - 1;

global OmegaM;  % IS THIS THE RIGHT OMEGA FOR THE GREEN FUNCTION?
global scalefactor;

% Determine the Green's function array for calculating phi
for l = 0:gridSize-1
    for m = 0:gridSize-1
        for n = 0:gridSize-1
            if (l==0)&&(m==0)&&(n==0)
                Green(1,1,1)=0;
            else
                Green(l+1,m+1,n+1)= ...
                    -(3*OmegaM/(8*scalefactor)) * ( sin(pi*l/gridSize)^2 + sin(pi*m/gridSize)^2 + sin(pi*n/gridSize)^2 )^(-1);
            end
        end
    end
end

% Solve the discretized Poisson Equation
% Need to solve in k-space where the equation becomes phik = Green*deltak
% Get the overdensity in k-space using FFT

DeltaK = fftn(DeltaX);

% Compute Phi in fourier space using the Green function
for l = 1:gridSize
    for m = 1:gridSize
        for n = 1:gridSize
            PhiK(l,m,n)= Green(l,m,n)*DeltaK(l,m,n);
        end
    end
end

% Now we FFT phik back to configuration space to get phi
Phi = real(ifftn(PhiK));
