function rocket_open_loop
clear all
close all

% Define system parameters
%
g = 9.81;                   % gravitational acceleration in m/s^2
mass = 5;                   % mass in (kg)
L = 10;                     % length of rocket (m)
thrust = mass*g*1.1;        % thrust of rocket (N)
I = (1/12)*mass*L^2;        % mass moment of inertia (kg m^2);

% Specify magnitude and angle of initial thrust vector
 
theta = 0.05;                  % initial angle (degrees) of thrust vector
thrust = mass*g*2.1;

% Specify initial velocity. Inital displacements are taken to be zero. 
z1_0 = 0;       % init x pos
z2_0 = 0;		% init x vel
z3_0 = L/2;    	% init y pos 				
z4_0 = 0;		% init y vel
z5_0 = 0;    	% init phi pos 				
z6_0 = 0;		% init phi dot

% put IC's into a single variable 
Z_0 = [z1_0, z2_0, z3_0, z4_0, z5_0, z6_0];

% Define simulation parameters
T_span = [0: 0.05: 30];  % time range for simulation with specified time step

% integration will stop at end of time span or when projectile hits ground
% (as defined by events)

options = odeset('Events', @event_stop);
[t, zout] = ode45(@projectile_fun, T_span, Z_0, options);

assignin('base','t',t);
assignin('base','zout',zout);
t(end)  %show end time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % The following section plots the rocket
figure(1)
for time=1:int16(length(t)/(20*max(t))):length(t)
    % Compute the tip and bottom of the rocket
    rocketTop = [zout(time, 1) + (L/2)*sin(-zout(time, 5)), zout(time, 3) + (L/2)*cos(-zout(time, 5))];
    rocketBot = [zout(time, 1) - (L/2)*sin(-zout(time, 5)), zout(time, 3) - (L/2)*cos(-zout(time, 5))];
    % Plot the ground
    plot([-10 10],[0 0],'k','LineWidth',2), hold on
    % plot the rocket
    scatter([rocketBot(1),rocketTop(1)], [rocketBot(2),rocketTop(2)],[],[1,0,0; 0,0,0]); 
    plot([rocketBot(1),rocketTop(1)], [rocketBot(2),rocketTop(2)],'k', 'LineWidth', 2); 
    
    % Some window sizing functions
    axis([-2 2 -2 120]); axis equal
    grid on
    set(gcf,'Position',[100 100 1000 400])
    drawnow limitrate, hold off
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following section is the ODE45 dzdt function. 

function dzdt = projectile_fun(T,Z)

% solve simultaneously for x and y,  and dx/dt and dy/dt
% Z(1)= z1=x, Z(2)= z2=dx/dt, Z(3)= z3=y, Z(4)= z4=dy/dt

dz1dt = Z(2); 
dz2dt = -(1/mass)*(thrust*sin(theta+Z(5)));
dz3dt = Z(4);
dz4dt = (1/mass)*(thrust*cos(theta+Z(5))) - g;
dz5dt = Z(6);
dz6dt = (-L*thrust*sin(theta))/(2*I);

dzdt = [dz1dt;dz2dt;dz3dt;dz4dt;dz5dt;dz6dt]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stops the simulation if the rocket hits the ground

function [eventvalue,stopthecalc,eventdirection] = event_stop(T,Z)
     
 % stop when Z(3)= z3= y = 0 (mass hits the ground in y-dir)
        eventvalue      =  Z(3);    %  ‘Events’ are detected when eventvalue=0
        stopthecalc     =  1;       %  Stop if event occurs
        eventdirection  = -1;       %  Detect only events with dydt<0
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end