%set initial X,Y,Z positions based on power spectrum

%random for now
% X = gridSize*rand(particleNum,1);
% Y = gridSize*rand(particleNum,1);   %pNx1 array of random numbers between 0 and gridSize
% Z = zeros(particleNum,1); Z = Z + gridSize/2;
% 
% X = [30,35];
% Y = [35,35];   
% Z = zeros(particleNum,1); Z = Z + gridSize/2;
% % 
% %evenly spaced initial grid:
N = 1;
for x = 0:gridSize/2 - 1
    for y = 0:gridSize/2 - 1
         for z = 0:gridSize/2 - 1
            X(N)=2*x + 0.2*randn;    %gaussian noise
            Y(N)=2*y + 0.2*randn;
            Z(N)=2*z + 0.2*randn;
            N = N+1;
         end
    end
end
% 
%     