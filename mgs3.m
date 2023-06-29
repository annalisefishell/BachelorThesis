function [A,R] = mgs3(A)
%
%   For a 3-tensor A, apply modified Gram-Schmidt to A(:,:,j) for each j.
%
%   Provided by Svetlana Dubinkina
%

[M,p,N] = size(A);
R = zeros(p,p,N);
for j=1:p
    R(j,j,:) = sqrt(sum(A(:,j,:).^2,1));
    A(:,j,:) = A(:,j,:)./repmat(R(j,j,:),[M,1]);
    R(j,j+1:p,:) = sum(repmat(A(:,j,:),[1,p-j,1]).*A(:,j+1:p,:),1);
    A(:,j+1:p,:) = A(:,j+1:p,:)-repmat(A(:,j,:),[1,p-j,1]).*repmat(R(j,j+1:p,:),[M,1,1]);
end