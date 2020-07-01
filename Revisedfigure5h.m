clear;clc
close all
% load('BalancedSmallGridBoundary.mat')
% figure;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx(101:end)/10^5/2.5,'r')
% load('NonbalancedSmallGridBoundary29Hzkh061.mat')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx(101:end)/10^5/2.5,'k')


load('BalancedSmallGridBoundary.mat')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx2(101:end)/10^5/2.5+3.5,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx3(101:end)/10^5/2.5+2.5,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx4(101:end)/10^5/2.5+1.3,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx5(101:end)/10^5/2.5+0.3,'r')
hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx6(101:end)/10^5/2.5,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx7(101:end)/10^5/2.5-0.75,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramV8(101:end)/10^5/2.5-1.2,'r')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx9(101:end)/10^5/2.5-2,'r')
% 
load('NonbalancedSmallGridBoundary29Hzkh061.mat')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx(101:end)/10^5/2.5+4,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx2(101:end)/10^5/2.5+3.5,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx3(101:end)/10^5/2.5+2.5,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx4(101:end)/10^5/2.5+1.3,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx5(101:end)/10^5/2.5+0.3,'k')
hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx6(101:end)/10^5/2.5,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx7(101:end)/10^5/2.5-0.75,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramV8(101:end)/10^5/2.5-1.2,'k')
% hold on;plot([101:6590+100]*0.5,2.4*4*10^8/400*seismogramVx9(101:end)/10^5/2.5-2,'k')

axis([101/2 3340 -0.8 0.44])
legend('Balanced SGFD Method','Non-Balanced SGFD method Vx')
xlabel('Travel time(ms)')
ylabel('Amp');
grid on
box on