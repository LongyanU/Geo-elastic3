
clear;clc
close all
load('balancedM7.mat')
coeff
figure;plot(seismogramVz5(107:end)/10^5/2.5,'r')
load('NonbalancedM7.mat')
hold on;plot(seismogramVz5(107:end)/10^5/2.5,'k')
coeff


% load('BalancedSmallGridBoundary.mat')
% hold on;plot(seismogramVz2(101:end)/10^5/2.5+1,'r')
% hold on;plot(seismogramVz3(101:end)/10^5/2.5+0.75,'r')
% hold on;plot(seismogramVz4(101:end)/10^5/2.5+0.35,'r')
% hold on;plot(seismogramVz5(101:end)/10^5/2.5,'r')
% hold on;plot(seismogramVz6(101:end)/10^5/2.5-0.25,'r')
% hold on;plot(seismogramVz7(101:end)/10^5/2.5-0.5,'r')
% hold on;plot(seismogramVz8(101:end)/10^5/2.5-0.75,'r')
% hold on;plot(seismogramVz9(101:end)/10^5/2.5-1,'r')
% 
% load('NonbalancedSmallGridBoundary29Hzkh061.mat')
% hold on;plot(seismogramVz2(101:end)/10^5/2.5+1,'k')
% hold on;plot(seismogramVz3(101:end)/10^5/2.5+0.75,'k')
% hold on;plot(seismogramVz4(101:end)/10^5/2.5+0.35,'k')
% hold on;plot(seismogramVz5(101:end)/10^5/2.5,'k')
% hold on;plot(seismogramVz6(101:end)/10^5/2.5-0.25,'k')
% hold on;plot(seismogramVz7(101:end)/10^5/2.5-0.5,'k')
% hold on;plot(seismogramVz8(101:end)/10^5/2.5-0.75,'k')
% hold on;plot(seismogramVz9(101:end)/10^5/2.5-1,'k')

axis([107 5300 -1.7 1])
legend('Balanced SGFD Method','Non-Balanced SGFD method Vz')
xlabel('Travel time(ms)')
ylabel('Amp');
grid on




