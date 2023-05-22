clc;
clear all;
close all;
%define and plot the input signal
a = 2;
t = 0:2*pi/50:2*pi;
x = a*sin(t);
l = length(x);
figure
plot(x,'r');
%define and plot the delta modulated signal
delta = 0.2;
hold on
xn = 0;
for i = 1:l
    if x(i)>=xn(i)
        d(i)=1;
        xn(i+1) = xn(i)+delta;
    else
        d(i)=0;
        xn(i+1)=xn(i)-delta;
    end
end
stairs(xn,'b')
hold on
legend('Original','Dm signal')
%recover the original signal(apply demodulation)
for i = 2:d
    if d(i)>xn(i)
        d(i)=0;
        xn(i+1) = xn(i)-delta;
    else
        d(i)=1;
        xn(i+1)= xn(i)+delta;
    end
end
figure
plot(x,'r');
hold on
stairs(xn,'b')
plot(xn,'c')
legend('Original','Dm signal','Recovered')
