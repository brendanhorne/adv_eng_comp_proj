clc, clear, close
syms x f e d c EI q L
v = f*x^5 + e*x^4 + d*x^3 + c*x^2 - (f*L^4 + e*L^3 + d*L^2 + c*L)*x;
d1v = diff(v,x);
d2v = diff(d1v,x);
Ut = 0.5 * EI * int((d2v)^2,x,0,L) + int(q*x/L*v,x,0,L);
dUtdf = diff(Ut,f);
dUtde = diff(Ut,e);
dUtdd = diff(Ut,d);
dUtdc = diff(Ut,c);
eqns =  [dUtdf dUtde dUtdd dUtdc];
S = solve(eqns,[f e d c]);
S.f;
S.e;
S.d;
S.c;
clear syms;
syms x EI q l;
f = S.f;
e = S.e;
d = S.d;
c = S.c;
v = subs(v)