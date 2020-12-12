% MPC control of thrust vectoring rocket
% We first define our rocket state space model and other parameters
clear

g = 9.81;                   % gravitational acceleration in m/s^2
mass = 5;                   % mass in (kg)
L = 10;                 % length of rocket (m)
thrust = mass*g*1.1;              % thrust of rocket (N)
I = (1/12)*mass*L^2;    % mass moment of inertia (kg m^2);
%
% Specify magnitude and angle of initial velocity vector 
theta = 0;     % initial angle (degrees) of thrust vector

% Fixed points
% thrust = mass*g
% z1,z2,z3,z4,z5,z6,theta = 0
% (-1/mass)*(thrust*cos(theta + z5))
F0 = mass*g;
angled = 0;
df2dz5 = (-1/mass)*(F0*cos(angled));
df4dz5 = (-1/mass)*(F0*sin(angled));

df2du1 = (-1/mass)*sin(angled);
df2du2 = (-1/mass)*(F0*cos(angled));
df4du1 = (1/mass)*cos(angled);
df4du2 = (-1/mass)*F0*sin(angled);
df6du1 = (-L*sin(angled))/(2*I);
df6du2 = (-L*cos(angled))/(2*I);

% Linearized A matrix. Is not a good representation because eigs of A are
% zero. 
Ac = [0 1 0 0 0 0; 
     0 0 0 0 df2dz5 0;
     0 0 0 1 0 0;
     0 0 0 0 df4dz5 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0];
 
% Linearized B matrix
Bc = [0 0;
     df2du1 df2du2;
     0 0;
     df4du1 df4du2;
     0 0;
     df6du1 df6du2];
 
% output C matrix. Assuming full state feedback for the sake of simplicity
Cc = eye(6);

% Feed through matrix D. Is zero because we assume that the current input
% cannot affect the current output. This is an implicity assumption in
% receding horizon control. 

Dc = zeros(6,2);
 
% Check how controllable the system is
% rank(ctrb(A,B)) % if the rank is 6, it is fully controllable

% The next step is to make the system discrete
% Covert the continuous model to a discrete model by sampling every 0.5
% seconds
Delta_t = 0.5;
[Ad, Bd, Cd, Dd] = c2dm(Ac, Bc, Cc, Dc, Delta_t);

% define control and prediction horizons
% define tracking signal. in this case, get to zero
Nc = 5;
Np = 40;
rs = [0,0,0,0,0,0];
% cost of actuation
rw = 0;

% Once we have the augemented state space model, we can begin building the MPC
% controller

[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mimompcgain(Ad, Bd, Cd, Nc, Np, rs);
% [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mimompcgain(0.8, 0.1, 1, 4, 10, 1);

% define initial state vector
x_0 = [1;1; 1;1; 1;1; 1;1; 1;1; 1;1];
% x_0_test = [0.1; 0.2];
deltaU = (Phi_Phi+rw*eye(size(Phi_Phi)))\(Phi_R - Phi_F*x_0);

% The next step is to find the closed loop stability


