clear, clc

%% Program values
dt = 0.01;
tspan = 500;
T = [0 tspan];
timeframe = 0:dt:tspan;

%% Lorenz '96 values
% Paremeters
d = 18;
F = 8;

% Initial condition
Z0 = zeros(d,1);
for i=1:d
    Z0(i) = sin(2*pi*((i-1)/d));
end

%% Lyapunov exponents + Estimations
Z = zeros(d,length(timeframe));
Z(:,1) = Z0 + dt * randn(size(Z0));

% Fundamental matrix
XN = zeros(d,d,length(timeframe));
XN(:,:,1) = eye(d); 
[Q, R] = mgs3(XN(:,:,1));

AQ = zeros(d,d,length(timeframe));

for i = 1:length(timeframe)-1
    % Estimates of Lorenz 96
    [z, AN] = Lorenz96(Z(:,i), Q(:,:,i), F, d);
    % [z, AN] = Lorenz96(Z(:,i), XN(:,:,i), F, d);
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
xlim([0,20000])

for i = 1:d
    lambdas = log(squeeze(R(i,i,:)))/dt;

    % Calculate LEs
    final_lambda(i) = mean(lambdas);

    % Plot LEs
    plot(lambdas)
    yline(final_lambda(i), 'LineWidth', 2)
end
