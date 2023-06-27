function [xdot,Av] = Lorenz63(x,v,sigma,beta,rho)
%
%   [xdot,Av] = fL63(x,v,sigma,beta,rho);
%
%   Computes the vector field xdot = f(x), and the Jacobian product
%   Av = A*v, where A = f'(x)
% 

xdot = [sigma*(x(2,:)-x(1,:));
    x(1,:).*(rho - x(3,:)) - x(2,:);
    x(1,:).*x(2,:) - beta*x(3,:)]; % Lorenz 63 function

    x = reshape(x,size(x,1),1,size(x,2));
    Av = [ sigma*(v(2,:,:)-v(1,:,:));
        v(1,:,:).*repmat(rho - x(3,:,:),[1,3,1]) - repmat(x(1,:,:),[1,3,1]).*v(3,:,:) - v(2,:,:);
        repmat(x(1,:,:),[1,3,1]).*v(2,:,:) + v(1,:,:).*repmat(x(2,:,:),[1,3,1]) - beta*v(3,:,:)]; % A*v
    % what is p, expected dimensions of x,v
    % p - num LE