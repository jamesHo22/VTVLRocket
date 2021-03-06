function [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np)
% input to the function is the state space model.
% Ap, Bp, Cp: state space model
% Nc, control horizon
% Np, prediction horizon

% Create the augmented state space model
[m1, n1] = size(Cp);
[n1, n_in] = size(Bp);
A_e = eye(n1+m1, n1+m1);
A_e(1:n1, 1:n1) = Ap;
A_e(n1+1:n1+m1, 1:n1) = Cp*Ap;
B_e = zeros(n1+m1, n_in);
B_e(1:n1, :) = Bp;
B_e(n1+1:n1+m1,:) = Cp*Bp;
C_e = zeros(m1, n1+m1);
C_e(:, n1+1:n1+m1) = eye(m1, m1);


% compute the F and Phi matricies
n=n1+m1;
h(1,:) = C_e;
F(1,:) = C_e*A_e;
for kk=2:Np
    h(kk, :) = h(kk-1, :)*A_e;
    F(kk, :) = F(kk-1, :)*A_e;
end

v = h*B_e;
Phi = zeros(Np, Nc);
Phi(:, 1) = v;

for i=2:Nc
    Phi(:, i) = [zeros(i-1, 1); v(1:Np-i+1, 1)]; % Toeplotz matrix
end

% Assumes r(ki) = 1. Set-point is equal to one
% Assumes x(ki) = [0.1 0.2]'

BarRs = ones(Np, 1)*1; % set-point is one. [1 1 ... 1]*r(ki)
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;
Phi_R = Phi'*BarRs;
% 
% % Calculate the optimal control signal given rw
% rw = 0;
% inv(Phi_Phi + rw*eye(4,4))*(Phi_R - Phi_F*[0.1 0.2]')





