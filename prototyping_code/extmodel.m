% Create a simple continous model
Ac = [0 1 0; 3 0 1; 0 1 0];
Bc = [1; 1; 3];
Cc = [0 1 0];
Dc = zeros(1,1);

% Covert the continuous model to a discrete model by sampling every 1
% second
Delta_t = 1;
[Ad, Bd, Cd, Dd] = c2dm(Ac, Bc, Cc, Dc, Delta_t);

% Create the augmented state space model
[m1, n1] = size(Cd);
[n1, n_in] = size(Bd);
A_e = eye(n1+m1, n1+m1);
A_e(1:n1, 1:n1) = Ad;
A_e(n1+1:n1+m1, 1:n1) = Cd*Ad;
B_e = zeros(n1+m1, n_in);
Be(1:n1, :) = Bd;
Be(n1+1:n1+m1, :) = Cd*Bd;
C_e = zeros(m1, n1+m1);
C_e(:, n1+1:n1+m1) = eye(m1, m1);