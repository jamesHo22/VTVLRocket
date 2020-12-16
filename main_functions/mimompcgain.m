function [F, Phi, BarRs, Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mimompcgain(Ad, Bd, Cd, Nc, Np, rs)
% Ad: discrete A matrix
% Bd: discrete B matrix
% Cd: discrete C matrix
% Nc: control horizon
% Np: prediction horizon
% rs: setpoint


% compute the F and Phi matricies
% The first step is to predict what happens to our state as it steps into
% the future Np time steps. Np is the length of the prediction. 
% Nc is the control horizon. It is the length of the control trajectory.
% Nc is less than or equal to Np. 

% function is confirmed against textbook example

[m1, n1] = size(Cd); % number of rows of C is number of outputs, columns is equal to number of states
[n1, n_in] = size(Bd); % rows: number of states, columns: number of inputs
A_e = eye(n1+m1, n1+m1); % make augmented A matrix 
A_e(1:n1, 1:n1) = Ad; % top left corner is A_d
A_e(n1+1:n1+m1, 1:n1) = Cd*Ad; % bottom left corner is Cd*Ad
B_e = zeros(n1+m1, n_in); % Create augmented B matrix
B_e(1:n1, :) = Bd; % top half is Bd
B_e(n1+1:n1+m1, :) = Cd*Bd; % bottom half is Cd*Bd
C_e = zeros(m1, n1+m1); % Create augmented C matrix
C_e(:, n1+1:n1+m1) = eye(m1, m1); % Back half of matrix is identity matrix

F(1:m1, :) = C_e*A_e;
h(1:m1, :) = C_e;

for i = m1+1:m1:Np*m1
    F(i:i+m1-1, :) = F(i-m1:i-1, :)*A_e;
    h(i:i+m1-1, :) = h(i-m1:i-1, :)*A_e;
end

v = h*B_e;
[r_v, c_v] = size(v);

% create Phi matrix, also known as the Toeplotz matrix
counter = 1;
Phi = zeros(r_v, c_v*(Nc));
Phi(:, 1:c_v) = v;
for i = c_v+1:c_v:(Nc)*c_v
    Phi(:, i:i+c_v-1) = [zeros(m1*counter, c_v); v(1:r_v - counter*m1, :)]; % Toeplotz matrix
    counter = counter + 1;
end

% Define the matricies needed to compute the optimal control inputs
BarRs = ones(r_v,m1);
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;
Phi_R = Phi'*BarRs;