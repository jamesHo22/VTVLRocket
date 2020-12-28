% MPC control of thrust vectoring rocket
% We first define our rocket state space model and other parameters
% TODOs:
% 1. simulate the nonlinear system and use the linear system to control it
% 2. get MPC to track a reference
% 3. constrained MPC

clear all
close all


g = 9.81;                   % gravitational acceleration in m/s^2
mass = 5;                  % mass in (kg)
L = 10;                     % length of rocket (m)
thrust = mass*g*1.1;        % thrust of rocket (N)
I = (1/12)*mass*L^2;        % mass moment of inertia (kg m^2);

% Specify magnitude and angle of initial velocity vector 
theta = 0;                  % initial angle (degrees) of thrust vector

% Linearize the system about the following fixed points
% thrust = mass*g
% z1,z2,z3,z4,z5,z6,theta = 0

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

% Linearized A matrix. eigs(Ac) are all zero, meaning that it is marginally
% stable
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
 
% output C matrix. Assuming full state feedback for the sake of simplicity.
% Will eventually use a kalman filter to estimate this from sensors
Cc = [0 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 1 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 1;];
  
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

% Once we have the augemented state space model, we can begin building the MPC
% controller

% SISO Test case for debugging.
% ac = [1 1; 0 1];
% bc = [0.5; 1];
% cc = [1 0];
% xm = [0;0]; % inital state vairable for the plant SISO
% y = 0; % inital output SISO
% [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mimompcgain(ac, bc, cc, Nc, Np, 1); % SISO

Nc = 10; % control horizon
Np = 60; % prediction horizon
rs = [0,0,50,0,0,0]; % Define setpoint signal. We want the rocket to go to this pose

% Here, we comupute the gain matricies for the augmented state space model
[F, Phi, BarRs, Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mimompcgain(Ad, Bd, Cd, Nc, Np, rs); % MIMO
[n, n_in] = size(B_e); % define some useful dimensions

xm = [50;0;100;0;0;0]; % inital state variable for the plant
y = xm; % inital output is equal to the initial state
Xf = zeros(n,1); % inital state feedback variable
t_end = 60;
N_sim = t_end/Delta_t; % define number of time steps to run
r = ones(N_sim,1)*rs; % setpoint we want the controller to get to
u = 0; % u(k-1) = 0 the initial control signal
rw = 9000; % Control input penalty

% A limitation of this model is that it is an embedded integrator. It does
% not work if there are more outputs than inputs. Only works for square
% matricies. 

% trying different DeltaU equation to see if I can get it to perform like
% an LQR controller

% Q = eye(360);
% R = eye(20);
% rp = ones(360,6).*rs;
% 
% DeltaU_test = (R + Phi'*Q*Phi)\Phi'*Q*(rp - F*Xf);

for kk=1:N_sim 
    DeltaU = (Phi_Phi + rw*eye(size(Phi_Phi)))\(Phi'*BarRs*0-Phi_F*Xf);
    % DeltaU = (R + Phi'*Q*Phi)\Phi'*Q*(0 - F*Xf); % the other DeltaU
    % equation. 

    deltau = DeltaU(1:n_in,1);
    u = u+deltau;
    % Save the control output and state 
    u1(:, kk) = u;
    y1(:, kk) = y;
    xm_old = xm;
    % compute state
    % xm = ac*xm+bc*u; %SISO
    % This is using the linearized dyanmics, which is not accurate. I will
    % need to discritize the nonlinear dyanmics and use that here to
    % calculate the next state. 
    xm = Ad*xm +Bd*u; %MIMO
    % y = cc*xm; %siso
    y = Cd*xm; %mimo
    Xf = [xm-xm_old;y];
end

figure(1) 
subplot(411)
plot(y1(1,:))
legend('X position')
grid on
title('State Variables vs Sampling Interval')
subplot(412)
plot(y1(3,:))
legend('Y position')
grid on
subplot(413)
plot(u1(1,:))
legend('thrust')
grid on
subplot(414)
plot(u1(2,:))
legend('angle')
xlabel('sampling instant')
grid on

% % Animation Plots
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
