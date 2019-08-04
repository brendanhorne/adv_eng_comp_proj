clc, clear, close
syms x L T t EI m;
v = (3*x^2*L-x^3)*cos(((2*pi)/T)*t);
d1v = diff(v,x)
d2v = diff(d1v,x)
d2v = subs(d2v,t,0)
Smax = int(0.5 * EI * d2v^2,x,0,L)
v = subs(v,x,L)
v = subs(v,t,0)
Kmax = 0.5 * m * (((2*pi)/T)*v)^2
eqn = Smax == Kmax
eqn = subs(eqn,[L EI m],[1 200 100])
T = solve(eqn,T)
% T = subs(T,[L EI m],[1 200 100])
result = double(T(1))