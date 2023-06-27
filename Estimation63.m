clear, clc

%% Program values
% dimensions
d = 3;

% time
dt = 0.01;
tspan = 500;
T = [0 tspan];
timeframe = 0:dt:tspan;

%% Lorenz 63 values
% initial condition
Z0 = [1 1 1];

% constants
rho = 28;
sigma = 10;
beta = 8/3;

%% Estimations
% Euler
E = zeros(length(timeframe), d);
E(1,:) = Z0 + 0.01 * randn(size(Z0));

% Runge-Kutta
RK = zeros(length(timeframe), d);
RK(1,:) = Z0 + 0.01 * randn(size(Z0));

for i = 1:length(timeframe)-1
    % Euler
    [z,Az] = Lorenz63(E(i,:)', eye(d), sigma, beta, rho);
    E(i+1,:) = E(i,:) + dt * z';

    % Runge-Kutta
    RK(i+1,:) = RK4(@(z)Lorenz63(RK(i,:)', eye(d), sigma, beta, rho),...
        dt, RK(i,:)');

end

%% Plot the euler estimated model values
figure;
plot3(E(:,1), E(:,2), E(:,3));
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz 63 - Euler Method');
grid on;

%% Plot the runge-kutta estimated model values
figure;
plot3(RK(:,1), RK(:,2), RK(:,3));
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz 63 - Runge-Kutta');
grid on;