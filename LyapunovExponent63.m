clear, clc

%% Program values
% dimensions
d = 3;

% time
dt = 0.01;
tspan = 350;
T = [0 tspan];
timeframe = 0:dt:tspan;

%% Lorenz 63 values
% initial condition
Z0 = ones(d,1);

% constants
rho = 28;
sigma = 10;
beta = 8/3;

%% Lyapunov exponents + Estimations
Z = zeros(d,length(timeframe));
Z(:,1) = Z0 + dt * randn(size(Z0));

% Fundamental matrix
XN = zeros(d,d,length(timeframe));
XN(:,:,1) = eye(d); % random initial values (eye(d))
[Q, R] = mgs3(XN(:,:,1));

AQ = zeros(d,d,length(timeframe));

for i = 1:length(timeframe)-1
    % Estimates of Lorenz 63
    [z, AN] = Lorenz63(Z(:,i), Q(:,:,i), sigma, beta, rho);
    % [z, AN] = Lorenz63(Z(:,i), XN(:,:,i), sigma, beta, rho);
    Z(:,i+1) = Z(:,i) + dt * z;

    % Estimates of Q_n+1 and R_n+1
    AQ(:,:,i+1) = Q(:,:,i) + dt * AN;
    [Q(:,:,i+1), R(:,:,i+1)] = mgs3(AQ(:,:,i+1));
    
    % Ill-conditioned estimates
    % XN(:,:,i+1) = XN(:,:,i) + dt * AN;
    % [Q(:,:,i+1), R(:,:,i+1)] = mgs3(XN(:,:,i+1));
end


%% Estimates and plots of final Lyapunov exponents
final_lambda = zeros(d,1);

figure;
hold on;
xlim([0 length(timeframe)])
title("Approximating the LES")
ylabel("Estimates of LEs")
xlabel("Time frame")

for i = 1:d
    lambdas = log(squeeze(R(i,i,:)))/dt;

    % Calculate LEs
    final_lambda(i) = mean(lambdas);

    % Plot LEs
    plot(lambdas)
    yline(final_lambda(i), 'LineWidth', 2)
end

