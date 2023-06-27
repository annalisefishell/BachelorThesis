clear, clc

%% Program values
F = 8;

% dimensions
d = 18;

% time
dt = 0.01;
tspan = 500;
T = [0 tspan];
timeframe = 0:dt:tspan;

%% Lorenz '96 initial conditions
Z0 = zeros(d,1);
for i=1:d
    Z0(i) = sin(2*pi*((i-1)/d));
end

%% Estimations
% Euler
E = zeros(length(timeframe), d);
E(1,:) = Z0 + 0.01 * randn(size(Z0));

% Runge-Kutta
RK = zeros(length(timeframe), d);
RK(1,:) = Z0 + 0.01 * randn(size(Z0));

for i = 1:length(timeframe)-1
    % Euler
    [z,~] = Lorenz96(E(i,:),ones(d),F,d);
    E(i+1,:) = E(i,:) + dt * z';

    % Runge-Kutta
    RK(i+1,:) = RK4(@(z)Lorenz96(RK(i,:)',ones(d),F,d), dt, RK(i,:)');
end

%% Plot the euler estimated model values
figure;
plot3(E(:,1), E(:,2), E(:,3));
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz 96 - Euler Method');
grid on;

%% Plot the runge-kutta estimated model values
figure;
plot3(RK(:,1), RK(:,2), RK(:,3));
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz 96 - Runge-Kutta');
grid on;