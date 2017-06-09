function i = simpsons(f,a,b,m)
%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////
dx = (b-a)/(3*m);
weights = [1,2*ones(1,m-1)+2*mod((1:m-1),2),1] ;
x = a:(b-a)/m:b;
i = dx*f(x)*weights';
end