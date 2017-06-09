%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////
clear all;close all;clc

 %% Part 1
A = [1 1;0 1; 1 1];
[Q,R] = gram_schmidt(A);


[Q1,R1] = qr(A);
%% Part 2
clear all
clc
syms x
x = (0:11);
y = [47.2,53.2,60.6,70,79.4,89.3,91.7,88.9,82.4,71,56.9,47.7]';

a1 = [1;1;1;1;1;1;1;1;1;1;1;1];
a2 = cos(2*pi*x/12)';
a3 = sin(2*pi*x/12)';

Ap2 = [a1,a2,a3];
Aq = Ap2\y;

[Qp2,Rp2] = gram_schmidt(Ap2);

d = Qp2'*y;  %Calculates a values
Sol = Rp2\d; %Calculates a values
a1 = Sol(1);
a2 = Sol(2);
a3 = Sol(3);
 
Y = a1 + a2*cos(2*pi*x/12) + a3*sin(2*pi*x/12);

Predict = a1 + a2*cos(2*pi*(3.5)/12) + a3*sin(2*pi*(3.5)/12); %Temperature prediction for x=3.2 [April 15th]
figure
plot(x,y); hold on
plot(x,Y,'xr')
plot (3.5,Predict,'*g','MarkerSize',10,'MarkerFaceColor','g')
string = sprintf('Prediction: %g%c',Predict,char(176));
xlabel('Month');
ylabel('Temperature');
title('Temperature w/ Prediction')
legend('Data','Estimate',string)
hold off

%% Part 3,

clc
clear all

syms x
f1 = 1/(1+x^2);          %Part B
f2 = exp((sin(6*pi*x))); %Part C

a = -1;
b = 1;
f = matlabFunction(f1);
h = [1/2;1/4;1/8;1/16;1/32;1/64];

for i = 1:length(h)
   m(i) = (b-a)/h(i) ;
   I_mid(i) = midpoint(f,a,b,m(i));
   I_trap(i) = trapezoid(f,a,b,m(i));
   I_simp(i) = simpsons(f,a,b,m(i));
end


d2f = matlabFunction(diff(f1,x,2));
d4f = matlabFunction(diff(f1,x,4));


k = max(abs(d2f(a)),abs(d2f(b))); %Trap and Midpoint
k_simp = max(abs(d4f(a)),abs(d4f(b)));

for i = 1:length(h)
    error_mid(i) = k *(b-a)^3/(24*m(i)^2);
    error_trap(i) = k *(b-a)^3/(12*m(i)^2);
    error_simp(i) = k_simp *(b-a)^5/(180*m(i)^4);
end
figure
loglog(h,error_mid);hold on
loglog(h,error_trap)
loglog(h,error_simp)
hold off
xlabel('h')
ylabel('error')
title('Error Estimates')
legend('Midpoint','Trapezoid','Simpsons')

%% PArt C
clear all
clc
syms x
f1 = exp((sin(6*pi*x)));

f = matlabFunction(f1);
a = -1;
b = 1;
h = [1/2;1/4;1/8;1/16;1/32;1/64];

for i = 1:length(h)
   m(i) = (b-a)/h(i) ;
   I_mid(i) = midpoint(f,a,b,m(i));
   I_trap(i) = trapezoid(f,a,b,m(i));
   I_simp(i) = simpsons(f,a,b,m(i));
end



d2f = matlabFunction(diff(f1,x,2));
d4f = matlabFunction(diff(f1,x,4));


k = max(abs(d2f(a)),abs(d2f(b)));
k_simp = max(abs(d4f(a)),abs(d4f(b)));

for i = 1:length(h)
    error_mid(i) = k *(b-a)^3/(24*m(i)^2);
    error_trap(i) = k *(b-a)^3/(12*m(i)^2);
    error_simp(i) = k_simp *(b-a)^5/(180*m(i)^4);
end
figure
loglog(h,error_mid);hold on
loglog(h,error_trap)
loglog(h,error_simp)
hold off
xlabel('h')
ylabel('error')
title('Error Estimates')
legend('Midpoint','Trapezoid','Simpsons')


%% Check Function
f = @(x) 1./x;

m = 2;
a = 1;
b = 2;
cMid = midpoint(f,a,b,m); %(Should Get ~0.68571)
cTrap = trapezoid(f,a,b,m); %(Should Get ~0..708333)
cSimp = simpsons(f,a,b,m); %(Should Get ~0..69444)
d2f = matlabFunction(diff(f,x,2));

k = max(abs(d2f(a)),abs(d2f(b)));
error_mid = k *(b-a)^3/(24*m^2); %(Should Get ~ 1/48)
error_trap = k *(b-a)^3/(12*m^2); %(Should Get ~ 1/24)
error_simp = k *(b-a)^5/(180*m^4);

%% Check Error

a = 0;
b = 2;
m = [2,4,8,16,32,64,128];

syms x
f1 = exp(x.^2);

f = matlabFunction(f1);

for i = 1:length(m)
   I_mid(i) = midpoint(f,a,b,m(i));
   I_trap(i) = trapezoid(f,a,b,m(i));
   I_simp(i) = simpsons(f,a,b,m(i));
end

d2f = matlabFunction(diff(f1,x,2));
d4f = matlabFunction(diff(f1,x,4));

k = max(abs(d2f(a)),abs(d2f(b)));
k_simp = max(abs(d4f(a)),abs(d4f(b)));

for i = 1:length(m)
    error_mid(i) = k *(b-a)^3/(24*m(i)^2);
    error_trap(i) = k *(b-a)^3/(12*m(i)^2);
    error_simp(i) = k_simp *(b-a)^5/(180*m(i)^4);
end
h = (b-a)./m;
figure
loglog(h,error_mid);hold on
loglog(h,error_trap)
loglog(h,error_simp)
hold off
xlabel('h')
ylabel('error')
title('Error Estimates')
legend('Midpoint','Trapezoid','Simpsons')


