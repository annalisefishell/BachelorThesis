function zout = RK4(fun, dt, yk)
% RUNGEKUTTA produces the approximations of the function using the RK4
% method for one step
%
%   The input is the function (fun) being solved, the time step between each
%   approximations (dt) and the guess from one time step previously (yk)
%
%   The output is the approximation at the time set (zout)
%
k1 = fun(yk);
k2 = fun(yk + (dt/2).*k1);
k3 = fun(yk + (dt/2).*k2);
k4 = fun(yk + dt.*k3);

zout = yk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);