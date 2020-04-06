
clear;clc
close all
load('BalancedSmallGridBoundary.mat')
figure;plot(2.4*4*10^8/400*seismogramVz(101:end)/10^5+1.5,'r')
load('NonbalancedSmallGridBoundary.mat')
hold on;plot(2.4*4*10^8/400*seismogramVz(101:end)/10^5+1.5,'k')
legend('Non-Balanced SGFD method Vz','Balanced SGFD Method')

load('BalancedSmallGridBoundary.mat')
hold on;plot(2.4*4*10^8/400*seismogramVz2(101:end)/10^5+1,'r')
hold on;plot(2.4*4*10^8/400*seismogramVz3(101:end)/10^5+0.75,'r')
hold on;plot(2.4*4*10^8/400*seismogramVz4(101:end)/10^5+0.35,'r')
hold on;plot(2.4*4*10^8/400*seismogramVz5(101:end)/10^5,'r')
hold on;plot(2.4*4*10^8/400*seismogramVz6(101:end)/10^5-0.25,'r')
hold on;plot(2.4*4*10^8/400*seismogramVz7(101:end)/10^5-0.5,'r')
hold on;plot(2.4*4*10^8/400*seismogramVz8(101:end)/10^5-0.75,'r')
hold on;plot(2.4*4*10^8/400*seismogramVz9(101:end)/10^5-1,'r')

load('NonbalancedSmallGridBoundary.mat')
hold on;plot(2.4*4*10^8/400*seismogramVz2(101:end)/10^5+1,'k')
hold on;plot(2.4*4*10^8/400*seismogramVz3(101:end)/10^5+0.75,'k')
hold on;plot(2.4*4*10^8/400*seismogramVz4(101:end)/10^5+0.35,'k')
hold on;plot(2.4*4*10^8/400*seismogramVz5(101:end)/10^5,'k')
hold on;plot(2.4*4*10^8/400*seismogramVz6(101:end)/10^5-0.25,'k')
hold on;plot(2.4*4*10^8/400*seismogramVz7(101:end)/10^5-0.5,'k')
hold on;plot(2.4*4*10^8/400*seismogramVz8(101:end)/10^5-0.75,'k')
hold on;plot(2.4*4*10^8/400*seismogramVz9(101:end)/10^5-1,'k')

axis([0 3300 -1.5 2.75])
% legend('Non-Balanced FD method','Analytical Method')
xlabel('Travel time(ms)')
ylabel('Amp');
grid on




