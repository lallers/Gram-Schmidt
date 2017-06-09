function i = trapezoid(f,a,b,m)
%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////
dx = (b-a)/(m-1); 
x = linspace(a,b,m);                
y = f(x);                        
sum_y = sum(y)-(1/2)*(y(1) + y(end)) ; 
i = sum_y * dx;
end

