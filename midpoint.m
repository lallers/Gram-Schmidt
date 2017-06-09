function i = midpoint(f,a,b,m)
%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////
dx = (b-a)/m;
k = 1:m;
x = a + (1/2)*(2*k-1)*dx;
y = f(x);
i = dx*sum(y);
end
