function [xdot, Av] = Lorenz96(x,v,F,N)
%LORENZ96 implements the Lorenz 96 model in either first or second order
%   With constant forcing F and dimension N

% calculate vector field
xdot = zeros(N,1);

for i=3:N-1
    xdot(i) = ((x(i + 1) - x(i - 2)) .* x(i - 1) - x(i)) + F;
end
% edge cases with indexing
xdot(1) = ((x(2) - x(N-1)) * x(N) - x(1)) + F;
xdot(2) = ((x(3) - x(N)) * x(1) - x(2)) + F;
xdot(N) = ((x(1) - x(N - 2)) * x(N - 1) - x(N)) + F;


% Calculate the jacobian times a matrix/vector
x = reshape(x,size(x,1),1,size(x,2));

diagonal = zeros(N,3);

% edge cases with indexing
diagonal(1,1) = x(N);

diagonal(1,2) = -x(N-1)+x(2); diagonal(2,2) = -x(N)+x(3); 
diagonal(N,2) = -x(N-2)+x(1); 

diagonal(1,3) = -x(N); 

for i = 2:N
    diagonal(i,1) = x(i-1);
    diagonal(i,3) = -x(i-1);
    if i>2 && i<N
        diagonal(i,2) = -x(i-2)+x(i+1);
    end
end

A = diag(-1*ones(N,1)) + diag(diagonal(1:N-1,1),1) + ...
    diag(diagonal(2:N,2),-1) + diag(diagonal(3:N,3),-2) + ...
    diag(diagonal(1:2,3),N-2);

A(N,1) = diagonal(N,1);
A(1,N) = diagonal(1,2);

Av = A*v;