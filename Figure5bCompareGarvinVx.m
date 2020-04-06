% clear;clc
% close all


% load('NonBalanced4March19.mat')
% load('Balanced4March192.mat')
% load('NonBalanced4March1922.mat')
% load('NonBalancedChangeXPosition.mat')
% load('NonBalanced700point35.mat')
% load('NonBalanced7.mat')
% load('NonBalanced900point45.mat')

Garvin2(0.5,(101-46)*10/1000)
load('GarvinResult.mat')
figure;plot(2.4*4*10^8/400*seismogramVx(101:end)/10^5+4,'r')
temp=conv(diff(src),u)/400;
hold on; plot(temp(101:end)/10^5+4,'b')
% axis([0 3300 -1*10^5 1.3*10^5 ])
legend('Non-Balanced FD method Vx','Analytical Method')
% xlabel('Travel time(ms)')
% ylabel('Amp');
% grid on


Garvin2(1,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVx2(101:end)/10^5+3.5,'r')
temp=conv(diff(src),u)/400;
hold on; plot(temp(101:end)/10^5+3.5,'b')
% axis([0 3300 -0.5*10^5 0.7*10^5 ])
% legend('Non-Balanced FD method Vz','Analytical Method')
% xlabel('Travel time(ms)')
% ylabel('Amp');
% grid on


Garvin2(2,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;plot(2.4*4*10^8/400*seismogramVx3(101:end)/10^5+2.5,'r')
temp=conv(diff(src),u)/400;
hold on; plot(temp(101:end)/10^5+2.5,'b')
% axis([0 3300 -0.4*10^5 0.6*10^5 ])
% legend('Non-Balanced FD method Vz','Analytical Method')
% xlabel('Travel time(ms)')
% ylabel('Amp');
% grid on



Garvin2(3,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVx4(101:end)/10^5+1.3,'r')
temp=conv(diff(src),u)/400;
hold on; plot(temp(101:end)/10^5+1.3,'b')
% % axis([0 3300 -0.3*10^5 0.45*10^5 ])
% % legend('Non-Balanced FD method Vz','Analytical Method')
% % xlabel('Travel time(ms)')
% % ylabel('Amp');
% % grid on
% 
Garvin2(4,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVx5(101:end)/10^5+0.3,'r')
temp=conv(diff(src),u)/400;
hold on; plot(temp(101:end)/10^5+0.3,'b')
% % axis([0 3300 -0.25*10^5 0.4*10^5 ])
% % legend('Non-Balanced FD method Vz','Analytical Method')
% % xlabel('Travel time(ms)')
% % ylabel('Amp');
% % grid on
% 
% 
Garvin2(4.5,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVx6(101:end)/10^5-0.2,'r')
temp=conv(diff(src),u)/400;
hold on; plot(temp(101:end)/10^5-0.2,'b')
% % axis([0 3300 -0.2*10^5 0.3*10^5 ])
% % legend('Non-Balanced FD method Vz','Analytical Method')
% % xlabel('Travel time(ms)')
% % ylabel('Amp');
% % grid on
% 
% %figure 7
Garvin2(2,(150-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVx7(101:end)/10^5-0.75,'r')
temp=conv(diff(src),u)/400;
hold on; plot(temp(101:end)/10^5-0.75,'b')
% % axis([0 3300 -0.5*10^5 0.52*10^5 ])
% % legend('Non-Balanced FD method Vz','Analytical Method')
% % xlabel('Travel time(ms)')
% % ylabel('Amp');
% % grid on
% 
% %figure 8
Garvin2(2,(200-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramV8(101:end)/10^5-1.2,'r')
temp=conv(diff(src),u)/400;
hold on; plot(temp(101:end)/10^5-1.2,'b')
% % axis([0 3300 -0.5*10^5 0.52*10^5 ])
% % legend('Non-Balanced FD method','Analytical Method')
% % xlabel('Travel time(ms)')
% % ylabel('Amp');
% % grid on
% 
% % figure 9
Garvin2(2,2.04)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVx9(101:end)/10^5-2,'r')
temp=conv(diff(src),u)/400;
hold on; plot(temp(101:end)/10^5-2,'b')
% axis([0 3300 -3 5])
axis([0 3300 -3.5 5.2])
% legend('Non-Balanced FD method','Analytical Method')
xlabel('Travel time(ms)')
ylabel('Amp');
grid on
