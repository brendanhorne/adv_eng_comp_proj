clc, clear, close
% Question 2
% Write a script to determine how many terms in the series exp(1/k), k =
% 1,2,3..., are required for the sum of the terms to exceed 2000. What is
% the sum for this number of terms.
s=0;
k=0;

while s < 2000
    k = k+1;
    s = s + exp(1/k);
    S(k) = exp(1/k);
end

k
s

x = 1 : 1: k;
figure
plot(x,S)
title('k vs exp(1/k)')
xlabel('k')
ylabel('exp(1/k)')