% Step by step  Page 46 in book
% 1. create the plant. This one is already discrete

Ap = [1 1; 0 1];
Bp = [0.5; 1];
Cp = [1 0];
Dp = 0;
Np = 20;
Nc = 4;

% The program calls the function mpcgain.m to generate the necessary gain 
% matrices and specifies the initial conditions for implementation of 
% receding horizon control. The initial state variable for the plant is 
% xm=0; and the initial state feedback variable is Xf=0; set-point signal 
% is specified and the number of simulation points is specified as 100. 

[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Ap, Bp, Cp, Nc, Np);
[n, n_in] = size(B_e);
xm = [0;0];
Xf = zeros(n,1);
N_sim = 100;
r = zeros(N_sim,1);
r((N_sim/2):end) = 1;
u = 0; % u(k-1) = 0 
y = 0;
rw = 10;

for kk=1:N_sim
    DeltaU = inv(Phi_Phi + rw*eye(Nc, Nc))*(Phi_R*r(kk)-Phi_F*Xf);
    deltau = DeltaU(1,1);
    u = u+deltau;
    % Save the control output and state 
    u1(kk) = u;
    y1(kk) = y;
    xm_old = xm;
    xm = Ap*xm+Bp*u;
    y = Cp*xm;
    Xf = [xm-xm_old;y];
end

k = 0:(N_sim-1);
figure 
subplot(211)
plot(k,y1)
xlabel('Sampling Instant')
legend('Output')
subplot(212)
plot(k,u1)
xlabel('Sampling Instant')
legend('Control')
