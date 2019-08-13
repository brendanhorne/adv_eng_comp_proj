clc, clear, close
syms x d c EI q L
v = d*x^3 + c*x^2 - (d*L^2 + c*L)*x
d1v = diff(v,x)
d2v = diff(d1v,x)
Ut = 0.5 * EI * int((d2v)^2,x,0,L) + int(q*x/L*v,x,0,L)
dUtdd = diff(Ut,d)
dUtdc = diff(Ut,c)
eqns = [dUtdd dUtdc];
S = solve(eqns,[d c]);
S.d
S.c
clear syms;
syms x EI q l;
d = S.d
c = S.c
b = -d*L^2 - c*L
v = subs(v)
