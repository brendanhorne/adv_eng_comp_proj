syms theta L k P
TPE = 1/2*k*(1/2*L*theta)^2 - P*(1/2*L*theta^2);
d= diff(TPE,theta)
