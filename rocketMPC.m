% MPC control of thrust vectoring rocket
% We first define our rocket state space model and other parameters
clear all
close all

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

% The next step is to make the system discrete
% Covert the continuous model to a discrete model by sampling every 0.5
% seconds
Delta_t = 0.5;
[Ad, Bd, Cd, Dd] = c2dm(Ac, Bc, Cc, Dc, Delta_t);

% Define some useful dimensions
% [m1, n1] = size(Cd);
% [n1, n_in] = size(Bd);

Nc = 10; % control horizon
Np = 60; % prediction horizon

% Define setpoint signal. We want the rocket to go to this pose
rs = [0,0,0,0,0,5];

% Once we have the augemented state space model, we can begin building the MPC
% controller

% Here, we comupute the gain matricies for the augmented state space model
% [Phi, Phi_Phi, Phi_F, A_e, B_e, C_e] = mimompcgain(Ad, Bd, Cd, Nc, Np);

% SISO Test case for debugging
ac = [1 1; 0 1];
bc = [0.5; 1];
cc = [1 0];

% using mimompcgain to calculate gains
% [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mimompcgain(ac, bc, cc, Nc, Np, 1); % SISO
[Phi, BarRs, Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mimompcgain(Ad, Bd, Cd, Nc, Np, rs); % MIMO
[n, n_in] = size(B_e);

xm = [20;0;50;0;0;0]; % inital state vairable for the plant MIMO
% xm = [0;0]; % inital state vairable for the plant SISO
Xf = zeros(n,1); % inital state feedback variable
t_end = 300;
N_sim = t_end/Delta_t; % define number of time steps to run


% r = zeros(N_sim,1); % setpoint
% r((N_sim/2):end) = 1; % the second half is ones
r = ones(N_sim,1)*rs; % setpoint
r_test = [50 0 100 0 0 0;
          50 0 100 0 0 0;
          50 0 100 0 0 0;
          50 0 100 0 0 0;
          50 0 100 0 0 0;
          50 0 100 0 0 0;];

u = 0; % u(k-1) = 0 the initial control signal
y = [20;0;50;0;0;0]; % inital output MIMO
% y = 0; % inital output SISO
rw = 0;

for kk=1:N_sim
    DeltaU = (Phi_Phi + rw*eye(size(Phi_Phi)))\(Phi'*BarRs*0-Phi_F*Xf);
    deltau = DeltaU(1:n_in,1);
    u = u+deltau;
    % Save the control output and state 
    u1(:, kk) = u;
    y1(:, kk) = y;
    xm_old = xm;
%     % compute state
%     xm = ac*xm+bc*u; %SISO
    xm = Ad*xm+Bd*u; %MIMO
%     y = cc*xm; %siso
    y = Cd*xm; %mimo
    Xf = [xm-xm_old;y];
end

figure(1) 
subplot(411)
plot(y1(1,:), 'o')
legend('X position')
subplot(412)
plot(y1(3,:), 'o')
legend('Y position')
subplot(413)
plot(u1(1,:), 'o')
legend('thrust')
subplot(414)
plot(u1(2,:), 'o')
legend('angle')
xlabel('sampling instant')

% % % Plots
% figure(2)
% zout = y1';
% for time=1:N_sim
%     rocketTop = [zout(time, 1) + (L/2)*sin(-zout(time, 5)), zout(time, 3) + (L/2)*cos(-zout(time, 5))];
%     rocketBot = [zout(time, 1) - (L/2)*sin(-zout(time, 5)), zout(time, 3) - (L/2)*cos(-zout(time, 5))];
% %     Plots the ground
%     plot([-10 10],[0 0],'k','LineWidth',2), hold on
% %     plots the rocket
%     scatter([rocketBot(1),rocketTop(1)], [rocketBot(2),rocketTop(2)],[],[1,0,0; 0,0,0]); 
%     plot([rocketBot(1),rocketTop(1)], [rocketBot(2),rocketTop(2)],'k', 'LineWidth', 2); 
%     
% %     Plot thrust vector
% %     u = K*zout(time,:)';
% %     u1 = u(1);
% %     u2 = u(2);
% %     dp = u1*[sin(-zout(time, 5)-u2), cos(-zout(time, 5)-u2)];
% %     quiver(rocketBot(1),rocketBot(2),dp(1),dp(2))
% %     set some window params
%     axis([-2 2 -2 120]); axis equal
%     grid on
%     set(gcf,'Position',[100 100 1000 400])
%     drawnow limitrate, hold off
%     
% end
