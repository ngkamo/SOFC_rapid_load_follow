function [f] = equi_ref(x,nCH4in,Tr,P)
A = nCH4in;
B = 2.5*A;

a = x(1)*A;
b = A*x(1)*(1-x(2));

Kp1 = 1.198e13*exp(-26830/(Tr));
Kp2 = 1.767e-2*exp(-4400/(Tr));


f = 1e3*[Kp1/P^2 - ((a-b)*(3*a+b)^3)/((A-a)*(B-a-b)*(A+B+2*a)^2);
     Kp2 - b*(3*a-b)/((a-b)*(B-a-b))];
