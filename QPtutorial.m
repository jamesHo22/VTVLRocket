% Quadratic programming test/tutorial
[x1, x2] = meshgrid(0:.1:10, 0:.1:10);

cost = 0.4*x1.^2 - 5*x1 + x2.^2 - 6*x2 + 50;

% plot cost function
figure(1)
contour(x1, x2, cost, 10, 'LineWidth', 2)
xlabel('x_{1}')
ylabel('x_{2}')
zlabel('Cost')

% plot constraint lines
hold on
plot([0;8], [2;10], 'k', 'LineWidth', 2)
plot([0;10], [8;5], 'k', 'LineWidth', 2)

H = 2*[0.4 0; 0 1];
f = [-5; -6];
A = [1 -1; -0.3 -1];
b = [-2; -8];
LB = [0;0];
UB = [10;10];

[X,J] = quadprog(H, f, A, b, [],[], LB, UB);
costVal = J + 50;

plot(X(1), X(2), 'r.', 'MarkerSize', 30)