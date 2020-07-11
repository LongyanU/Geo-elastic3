clear;clc
% close all
load('balancedM7.mat')
figure;plot(seismogramVx(101:end)/10^5/2.5+1.1,'r')
load('NonbalancedM7.mat')
hold on;plot(seismogramVx(101:end)/10^5/2.5+1.1,'k')


load('balancedM7.mat')
% hold on;plot(seismogramVx2(101:end)/10^5/2.5+3.5,'r')
hold on;plot(seismogramVx3(101:end)/10^5/2.5+.5,'r')
% hold on;plot(seismogramVx4(101:end)/10^5/2.5+1.3,'r')
hold on;plot(seismogramVx5(101:end)/10^5/2.5+0,'r')
% hold on;plot(seismogramVx6(101:end)/10^5/2.5-0.2,'r')
hold on;plot(seismogramVx7(101:end)/10^5/2.5-1,'r')
% hold on;plot(seismogramV8(101:end)/10^5/2.5-1.2,'r')
hold on;plot(seismogramVx9(101:end)/10^5/2.5-2,'r')

load('NonbalancedM7.mat')
% hold on;plot(seismogramVx(101:end)/10^5/2.5+4,'k')
% hold on;plot(seismogramVx2(101:end)/10^5/2.5+3.5,'k')
hold on;plot(seismogramVx3(101:end)/10^5/2.5+.5,'k')
% hold on;plot(seismogramVx4(101:end)/10^5/2.5+1.3,'k')
hold on;plot(seismogramVx5(101:end)/10^5/2.5+0,'k')
% hold on;plot(seismogramVx6(101:end)/10^5/2.5-0.2,'k')
hold on;plot(seismogramVx7(101:end)/10^5/2.5-1,'k')
% hold on;plot(seismogramV8(101:end)/10^5/2.5-1.2,'k')
hold on;plot(seismogramVx9(101:end)/10^5/2.5-2,'k')

% axis([101/2 3340 -3 5])
axis([107 5300 -4 5])
legend('Balanced SGFD Method','Non-Balanced SGFD method Vx')
xlabel('Travel time(ms)')
ylabel('Amp');
grid on