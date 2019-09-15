clc, clear, close
% Question 1
% Write a MATLAB Function to take a non-negative integer n, as the only
% input argument and calculates n-factorial denoted by n!.
x=1;
n = input("Calculate n! n = ");
if (n<0)
    disp("Please enter an integer greater than 0")
elseif mod(n,1) > 0
    disp("Please enter an integer")
else
    for i = 1:n
        x = x*i;
    end
end
x
