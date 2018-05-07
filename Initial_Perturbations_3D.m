
%{
Inputs:
cells/side = gridSize
objects/side = [obj_side]
initial density = [f(particleMass,gridSize)]
P(k) (or just amplitude, for linear)
%}

%{
Internal:
    Simulation Constants
average spaceing = [avg_spacing]
maximum k for number of objects/side = [max_k]
    Intermediate
X_dev = zeros(obj_side, obj_side, obj_side);
Y_dev = zeros(obj_side, obj_side, obj_side);
Z_dev = zeros(obj_side, obj_side, obj_side);
Px_0 = zeros(obj_side, obj_side, obj_side);
Py_0 = zeros(obj_side, obj_side, obj_side);
Pz_0 = zeros(obj_side, obj_side, obj_side);

    Counters
l_mode
m_mode
n_mode
%}

%{
Outputs:

Total Particle Number = particleNum
X
Y
Z
Px
Py
Pz
%}

% Declare Variables 
%{
mass per particle = particleMass
initial density = ?
P(k) = Amp = ??
%}
global gridSize;

%Sim Conditions:
% a_init = 0.1;
% gridSize = 64;
obj_side = 32;

%Outputs:
particleNum = obj_side ^ 3;
X = zeros(particleNum, 1);
Y = zeros(particleNum, 1);
Z = zeros(particleNum, 1); 
Px = zeros(particleNum, 1);
Py = zeros(particleNum, 1);
Pz = zeros(particleNum, 1);

%Internal Variables:
X_int = zeros(obj_side, obj_side, obj_side);
Y_int = zeros(obj_side, obj_side, obj_side);
Z_int = zeros(obj_side, obj_side, obj_side);
X_dev = zeros(obj_side, obj_side, obj_side);
Y_dev = zeros(obj_side, obj_side, obj_side);
Z_dev = zeros(obj_side, obj_side, obj_side);
Px_0 = zeros(obj_side, obj_side, obj_side);
Py_0 = zeros(obj_side, obj_side, obj_side);
Pz_0 = zeros(obj_side, obj_side, obj_side);

% Calculate Internal Constants

avg_spacing = gridSize / obj_side;
k_fundamental = (2 * pi) / gridSize;
max_n_k = floor(obj_side / 2);
%max_n_k = 4;


% Set up Initial Grid
disp("setting up initial conditions...")

for k_pos = 1:obj_side
    for j_pos = 1:obj_side
        for i_pos = 1:obj_side
            X_int(i_pos,j_pos,k_pos) = avg_spacing * (i_pos - 1/2);
            Y_int(i_pos,j_pos,k_pos) = avg_spacing * (j_pos - 1/2);
            Z_int(i_pos,j_pos,k_pos) = avg_spacing * (k_pos - 1/2);            
        end
    end
end

% Generate Random Waves

%V1 - Both Random (Noise about expected P(k))
a_rand = normrnd(0,1,[max_n_k,max_n_k,max_n_k]);
b_rand = normrnd(0,1,[max_n_k,max_n_k,max_n_k]);

%V2 - Exact P(k) Match
%a_rand = 2*rand - 1
%b_rand = (1 - a_rand ^ 2) ^ (1/2)

% Calculate Perturbation Grid

amplitude = 0.03;

for k_pos = 1:obj_side
    for j_pos = 1:obj_side
        for i_pos = 1:obj_side
            for n_mode = 1:max_n_k
                for m_mode = 1:max_n_k
                    for l_mode = 1:max_n_k
                        k_mag = (l_mode^2 + m_mode^2 + n_mode^2)^.5 * k_fundamental;
                        X_dev(i_pos,j_pos,k_pos) = X_dev(i_pos,j_pos,k_pos) + real(amplitude*(1i*(k_fundamental*l_mode)*(a_rand(l_mode,m_mode,n_mode)-b_rand(l_mode,m_mode,n_mode)*1i)/2)*(exp(1i*(k_fundamental*l_mode)*X_int(i_pos,j_pos,k_pos))));
                        Y_dev(i_pos,j_pos,k_pos) = Y_dev(i_pos,j_pos,k_pos) + real(amplitude*(1i*(k_fundamental*m_mode)*(a_rand(l_mode,m_mode,n_mode)-b_rand(l_mode,m_mode,n_mode)*1i)/2)*(exp(1i*(k_fundamental*m_mode)*Y_int(i_pos,j_pos,k_pos))));
                        Z_dev(i_pos,j_pos,k_pos) = Z_dev(i_pos,j_pos,k_pos) + real(amplitude*(1i*(k_fundamental*n_mode)*(a_rand(l_mode,m_mode,n_mode)-b_rand(l_mode,m_mode,n_mode)*1i)/2)*(exp(1i*(k_fundamental*n_mode)*Z_int(i_pos,j_pos,k_pos))));
                    end
                end
            end
        end
    end
end

% Calculate [Final] Initial Positions

for k_pos = 1:obj_side
    for j_pos = 1:obj_side
        for i_pos = 1:obj_side
            object_num = i_pos + obj_side * (j_pos - 1) + obj_side^2 * (k_pos - 1);
            
            X(object_num) = X_int(i_pos,j_pos,k_pos) + X_dev(i_pos,j_pos,k_pos);
            Y(object_num) = Y_int(i_pos,j_pos,k_pos) + Y_dev(i_pos,j_pos,k_pos);
            Z(object_num) = Z_int(i_pos,j_pos,k_pos) + Z_dev(i_pos,j_pos,k_pos);
            
            Px(object_num) = scalefactor * X_dev(i_pos,j_pos,k_pos);
            Py(object_num) = scalefactor * Y_dev(i_pos,j_pos,k_pos);
            Pz(object_num) = scalefactor * Z_dev(i_pos,j_pos,k_pos);    
        end
    end
end

% Display Particle Positions

%scatter(X_int(:,1,1), X_dev(:,1,1))
%scatter(X,Y)
%scatter3(X,Y,Z,'.');