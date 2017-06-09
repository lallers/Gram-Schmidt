function [Q,r] = gram_schmidt(A)
%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////

[n,m] = size(A);
Q = zeros(n,m);
r = zeros(m,m);

q(:,1) = A(:,1)/norm(A(:,1));
r(1,1) = norm(A(:,1));

for j = 2:m
    z = A(:,j);
    for i = 1:j-1
        r(i,j) = A(:,i)'*A(:,j);
        z = z - r(i,j)*q(:,i);
    end
    r(j,j) = norm(z);
    q(:,j) = z/r(j,j);
end
Q = q;
R = r;
end
