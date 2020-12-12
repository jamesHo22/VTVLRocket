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

% We then need to create the augmented state space model
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

% Once we have the augemented state space model, we can begin building the MPC
% controller

% The first step is to predict what happens to our state as it steps into
% the future Np time steps. Np is the length of the prediction. 
% Nc is the control horizon. It is the length of the control trajectory.
% Nc is less than or equal to Np. 

% Next, we will write our prediction equation that we will later optimize
% Lets start with just Np = Nc = 10. 

% compute the F and Phi matricies
Nc = 2;
Np = 10;

F(1:m1, :) = C_e*A_e;
h(1:m1, :) = C_e;

for i = m1+1:m1:Np*m1
    F(i:i+m1-1, :) = F(i-m1:i-1, :)*A_e;
    h(i:i+m1-1, :) = h(i-m1:i-1, :)*A_e;
end

v = h*B_e;
[r_v, c_v] = size(v);

% create Phi matrix
counter = 1;
Phi = zeros(r_v, c_v*(Nc));
Phi(:, 1:c_v) = v;
for i = c_v+1:c_v:(Nc)*c_v
    Phi(:, i:i+c_v-1) = [zeros(m1*counter, c_v); v(1:r_v - counter*m1, :)]; % Toeplotz matrix
    counter = counter + 1;
end

% optimizing (look at equation 1.16)
% Define setpoint signal
BarRs = ones(Np,m1);
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;
Phi_R = Phi'*BarRs;











