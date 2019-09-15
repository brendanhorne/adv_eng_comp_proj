clc, clear, close
% Using numerical integration with step of  ?x=0.02 find the integration
% above
syms x;
eqn = 2 + sqrt(x) / (x.+1);
x1 = 0:0.02:2;
T = 0;
for i = x1(1): length(x1)
    T = T + subs(eqn,x,i)*0.02;
end    
vpa(T)

B = traps(x1,eqn)