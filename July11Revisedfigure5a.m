
clear;clc
close all
load('balancedM4.mat')
figure;plot( seismogramVz(107:end)/10^5/2.5+2,'r')
load('NonbalancedM4.mat')
hold on;plot( seismogramVz(107:end)/10^5/2.5+2,'k')


load('balancedM4.mat')
% hold on;plot( seismogramVz2(107:end)/10^5/2.5+1.5,'r')
hold on;plot( seismogramVz3(107:end)/10^5/2.5+0.75,'r')
% hold on;plot( seismogramVz4(107:end)/10^5/2.5+0.9,'r')
hold on;plot( seismogramVz5(107:end)/10^5/2.5,'r')
% hold on;plot( seismogramVz6(107:end)/10^5/2.5-0.5,'r')
hold on;plot( seismogramVz7(107:end)/10^5/2.5-0.75,'r')
% hold on;plot( seismogramVz8(107:end)/10^5/2.5-0.9,'r')
hold on;plot( seismogramVz9(107:end)/10^5/2.5-1.3,'r')

load('NonbalancedM4.mat')
% hold on;plot( seismogramVz2(107:end)/10^5/2.5+1.5,'k')
hold on;plot( seismogramVz3(107:end)/10^5/2.5+0.75,'k')
% hold on;plot( seismogramVz4(107:end)/10^5/2.5+0.35,'k')
hold on;plot( seismogramVz5(107:end)/10^5/2.5,'k')
% hold on;plot( seismogramVz6(107:end)/10^5/2.5-0.5,'k')
hold on;plot( seismogramVz7(107:end)/10^5/2.5-0.75,'k')
% hold on;plot( seismogramVz8(107:end)/10^5/2.5-0.9,'k')
hold on;plot( seismogramVz9(107:end)/10^5/2.5-1.3,'k')

axis([107 5300 -4 5])
legend('Balanced SGFD Method','Non-Balanced SGFD method Vz')
xlabel('Travel time(ms)')
ylabel('Amp');
grid on




