

clc, clear, close
syms L landa;
L = 1;
K = [(288/5)*L^3 36*L^4 16*L^3; 36*L^4 24*L^3 12*L^2; 16*L^3 12*L^2 8*L];
KG = [(18/7)*L^7 2*L^6 (6/5)*L^5; 2*L^6 (8/5)*L^5 L^4; (6/5)*L^5 L^4 (2/3)*L^3];
Y = K / KG;
P = eig(Y)

