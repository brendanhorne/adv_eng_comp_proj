clc, clear, close
syms theta L k P lambda
TPE = 1/2*k*(1/2*L*theta)^2 - P*(1/2*L*theta^2)==0

TPE = TPE / theta^2
TPE = expand(TPE)
TPE = TPE / k
TPE = expand(TPE)
TPE = subs(TPE, P/k,lambda)
A = lhs(TPE)
A = eig(A)
lambda = solve(A,lambda)
clear syms
syms L k P
eqn = P/k == lambda
P = solve(eqn,P)