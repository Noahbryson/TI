close all
figure(1)
fs = 100000;
dt = 1/fs;
t = 0:dt:1;
base = 10000;
deltaF = 5;
% x1 = square(t*5000*2*pi);
x1 = sin(t*base*2*pi);
x2 = sin(t*(base+deltaF)*2*pi);
y = x1+x2;

env = hilbert(y);
% env2 = -sqrt(real(env).^2 + imag(env).^2);
% env2 = y - mean(y);
env2 = getLowPassData(abs(y'),1000,2,fs);
subplot(2,1,1)
plot(t,y);
hold on
plot(t,x1);
ylim([-3 3])
subplot(2,1,2)
plot(t,abs(env));
hold on
plot(t,env2)