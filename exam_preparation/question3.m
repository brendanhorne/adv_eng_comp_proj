clc, clear, close
% Draw 4 graphs of sin(x), sin2(x), sin3(x) and sin4(x) in one figure using
% 4 sub-plots. Use an interval of 0<x<2?

x = 0: pi/100: 2*pi;
figure

subplot(2,2,1);
plot(x,sin(x));
title('x vs sin(x)')
xlabel('x')
ylabel('sin(x)')

subplot(2,2,2);
plot(x,sin(x).^2);
title('x vs sin(x)^2')
xlabel('x')
ylabel('sin(x)^2')

subplot(2,2,3);
plot(x,sin(x).^3);
title('x vs sin(x)^3')
xlabel('x')
ylabel('sin(x)^3')

subplot(2,2,4);
plot(x,sin(x).^4);
title('x vs sin(x)^4')
xlabel('x')
ylabel('sin(x)^4')

figure
plot(x,sin(x),x,sin(x).^2,x,sin(x).^3,x,sin(x).^4);
legend('sin(x)','sin(x)^2','sin(x)^3','sin(x)^4');