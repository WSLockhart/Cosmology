% Script finding the potential phi By EED

% Calculating average density and overdensity
rhoav = mean(rho(:));
deltax = rho - 1;

% Start solving the discretized Poisson Equation
% Need to solve in k-space where the equation becomes phik = Gr*deltak
% Get the overdensity in k-space using FFT

deltak = fftn(deltax);

% Finding phi in fourier space
for l = 1:Ng
    for m = 1:Ng
        for n = 1:Ng
            phik(m,l,n)= Gr(l,m,n)*deltak(l,m,n);
        end
    end
end

% Now we FFT phik back to configuration space to get phi
phi = circshift(circshift(circshift(...  Shifting is to align the potential field to the mass positions
    real( ...   potential must be real
    ifftn(phik) ...
    ),-1,1),-1,2),-1,3);

%plotting bullshit
hold off
contourf(phi(:,:,1))
hold on
scatter3(x,y,z,'filled')
axis([0 Ng 0 Ng 0 Ng])
grid on