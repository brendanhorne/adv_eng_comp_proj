syms theta L k P
TPE = 1/2*k*(1/2*L*theta)^2 - P*(1/2*L*theta^2);
dTPEdtheta = diff(TPE,theta);
P = solve(dTPEdtheta,P);
P