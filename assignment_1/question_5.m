clc, clear, close
syms theta L k P lambda
TPE = 1/2*k*(1/2*L*theta)^2 - P*(1/2*L*theta^2)==0
TPE = TPE / theta^2
TPE = expand(TPE)
TPE = TPE / k
TPE = expand(TPE)
TPE = subs(TPE, (P*L)/(2*k),lambda*L/2)
A = lhs(TPE)
A = eig(A)
A = subs(A,lambda, P/k)
P = solve(A,P)