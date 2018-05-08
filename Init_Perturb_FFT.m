clearvars
close all;

%Sim Conditions:
gridSize = 128*2;
obj_side = 64*2;

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
X_c_k = zeros(obj_side + 1, obj_side + 1, obj_side + 1);
Y_c_k = zeros(obj_side + 1, obj_side + 1, obj_side + 1);
Z_c_k = zeros(obj_side + 1, obj_side + 1, obj_side + 1);
X_dev = zeros(obj_side, obj_side, obj_side);
Y_dev = zeros(obj_side, obj_side, obj_side);
Z_dev = zeros(obj_side, obj_side, obj_side);
Px_0 = zeros(obj_side, obj_side, obj_side);
Py_0 = zeros(obj_side, obj_side, obj_side);
Pz_0 = zeros(obj_side, obj_side, obj_side);


% Calculate Internal Constants

avg_spacing = gridSize / obj_side;
k_fund = (2 * pi) / gridSize;
max_n_k = floor(obj_side / 2);
%max_n_k = 32;


% Set up Initial Grid

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

a_rand = normrnd(0,1,[(2*max_n_k+1),(2*max_n_k+1),(2*max_n_k+1)]);
b_rand = normrnd(0,1,[(2*max_n_k+1),(2*max_n_k+1),(2*max_n_k+1)]);

Power_Spec = 1;
c_k(max_n_k+1,max_n_k+1,max_n_k+1) = 0;

for n = -max_n_k:max_n_k
    for m = -max_n_k:max_n_k
        for l = -max_n_k:max_n_k                       
            k_mag = (l^2 + m^2 + n^2)^.5 * k_fund;
   
            if k_mag ~= 0    
                c_k(l+max_n_k+1,m+max_n_k+1,n+max_n_k+1) = (k_mag)^(-2) * Power_Spec^(1/2) * (a_rand(l+max_n_k+1,m+max_n_k+1,n+max_n_k+1) - 1i * b_rand(l+max_n_k+1,m+max_n_k+1,n+max_n_k+1))/2;
            end

            X_c_k(l+max_n_k+1,m+max_n_k+1,n+max_n_k+1) = 1i * k_fund * l * c_k(l+max_n_k+1,m+max_n_k+1,n+max_n_k+1);
            Y_c_k(l+max_n_k+1,m+max_n_k+1,n+max_n_k+1) = 1i * k_fund * m * c_k(l+max_n_k+1,m+max_n_k+1,n+max_n_k+1);
            Z_c_k(l+max_n_k+1,m+max_n_k+1,n+max_n_k+1) = 1i * k_fund * n * c_k(l+max_n_k+1,m+max_n_k+1,n+max_n_k+1);
        end
    end
end

% Calculate Perturbation Grid

amplitude = 25;

X_dev = real(ifft(X_c_k,obj_side));
Y_dev = real(ifft(Y_c_k,obj_side));
Z_dev = real(ifft(Z_c_k,obj_side));

% Calculate [Final] Initial Positions

p_prefactor = 1;

for k_pos = 1:obj_side
    for j_pos = 1:obj_side
        for i_pos = 1:obj_side
            object_num = i_pos + obj_side * (j_pos - 1) + obj_side^2 * (k_pos - 1);
            
            X(object_num) = X_int(i_pos,j_pos,k_pos) + amplitude * X_dev(i_pos,j_pos,k_pos);
            Y(object_num) = Y_int(i_pos,j_pos,k_pos) + amplitude * Y_dev(i_pos,j_pos,k_pos);
            Z(object_num) = Z_int(i_pos,j_pos,k_pos) + amplitude * Z_dev(i_pos,j_pos,k_pos);  
            
            Px(object_num) = p_prefactor * X_dev(i_pos,j_pos,k_pos);
            Py(object_num) = p_prefactor * Y_dev(i_pos,j_pos,k_pos);
            Pz(object_num) = p_prefactor * Z_dev(i_pos,j_pos,k_pos);
        end
    end
end

% Display Particle Positions

scatter(X_int(:,1,1), amplitude * X_dev(:,1,1));
%scatter(X,Y);
%scatter3(X,Y,Z, '.');