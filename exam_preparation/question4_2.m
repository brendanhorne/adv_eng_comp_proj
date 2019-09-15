clc, clear, close
% Using the symbolic math find the exact value of the integration below: (10 points)
syms x;
eqn = 2 + sqrt(x) / (x+1);
A = int(eqn,0,2)
vpa(A)