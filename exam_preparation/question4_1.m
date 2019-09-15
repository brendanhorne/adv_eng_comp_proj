clc, clear, close
% Find the real roots of the polynomial x5-2x4+x2-x+1=0
syms x;
eqn = x^5 -2*x^4 + x^2 -x +1==0;
S = solve(eqn,x,'MaxDegree',5,'Real',true)
vpa(S)



