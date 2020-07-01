
clear;clc
close all
% load('BalancedSmallGridBoundary.mat')
% coeff
% figure;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz(101:end)/10^5/2.5,'r')
% load('NonbalancedSmallGridBoundary29Hzkh061.mat')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz(101:end)/10^5/2.5,'k')
% coeff


load('BalancedSmallGridBoundary.mat')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz2(101:end)/10^5/2.5+1,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz3(101:end)/10^5/2.5+0.75,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz4(101:end)/10^5/2.5+0.35,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz5(101:end)/10^5/2.5,'r')
hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz6(101:end)/10^5/2.5,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz7(101:end)/10^5/2.5-0.5,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz8(101:end)/10^5/2.5-0.75,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz9(101:end)/10^5/2.5-1,'r')
% 
load('NonbalancedSmallGridBoundary29Hzkh061.mat')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz2(101:end)/10^5/2.5+1,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz3(101:end)/10^5/2.5+0.75,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz4(101:end)/10^5/2.5+0.35,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz5(101:end)/10^5/2.5,'k')
hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz6(101:end)/10^5/2.5,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz7(101:end)/10^5/2.5-0.5,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz8(101:end)/10^5/2.5-0.75,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVz9(101:end)/10^5/2.5-1,'k')

axis([101/2 3300 -0.2 0.4])
legend('Balanced SGFD Method','Non-Balanced SGFD method Vz')
xlabel('Travel time(ms)')
ylabel('Amp');
grid on
box on




