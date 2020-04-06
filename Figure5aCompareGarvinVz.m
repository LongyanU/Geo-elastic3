
% load('March30Morning9.mat')
% load('NonBalanced7Vp800rho04txx0.mat')
Garvin2(0.5,(101-46)*10/1000)
load('GarvinResult.mat')
figure;plot(2.4*4*10^8/400*seismogramVz(101:end)/10^5+1.5,'r')
% figure;plot(2.4*4*10^8/400*seismogramVz(101:end)/10^5,'r')
temp=conv(diff(src),w)/400;
hold on; plot(temp(101:end)/10^5+1.5,'b')
% hold on; plot(temp(96:end)/10^5,'b')
% axis([0 3300 -1 1.3])
legend('Non-Balanced FD method Vz','Analytical Method')
xlabel('Travel time(ms)')
ylabel('Amp');
grid on


Garvin2(1,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;plot(2.4*4*10^8/400*seismogramVz2(101:end)/10^5+1,'r')
temp=conv(diff(src),w)/400;
hold on; plot(temp(101:end)/10^5+1,'b')
% % axis([0 3300 -0.5*10^5 0.7*10^5 ])
% % legend('Non-Balanced FD method Vz','Analytical Method')
% % xlabel('Travel time(ms)')
% % ylabel('Amp');
% grid on


Garvin2(2,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVz3(101:end)/10^5+0.75,'r')
temp=conv(diff(src),w)/400;
hold on; plot(temp(101:end)/10^5+0.75,'b')
% axis([0 3300 -0.4*10^5 0.6*10^5 ])
% legend('Non-Balanced FD method Vz','Analytical Method')
% xlabel('Travel time(ms)')
% ylabel('Amp');
% grid on



Garvin2(3,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVz4(101:end)/10^5+0.35,'r')
temp=conv(diff(src),w)/400;
hold on; plot(temp(101:end)/10^5+0.35,'b')
% axis([0 3300 -0.3*10^5 0.45*10^5 ])
% legend('Non-Balanced FD method Vz','Analytical Method')
% xlabel('Travel time(ms)')
% ylabel('Amp');
% grid on

Garvin2(4,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVz5(101:end)/10^5,'r')
temp=conv(diff(src),w)/400;
hold on; plot(temp(101:end)/10^5,'b')
% axis([0 3300 -0.25*10^5 0.4*10^5 ])
% legend('Non-Balanced FD method Vz','Analytical Method')
% xlabel('Travel time(ms)')
% ylabel('Amp');
% grid on
% 

Garvin2(4.5,(101-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVz6(101:end)/10^5-0.25,'r')
temp=conv(diff(src),w)/400;
hold on; plot(temp(101:end)/10^5-0.25,'b')
% % axis([0 3300 -0.2*10^5 0.3*10^5 ])
% % legend('Non-Balanced FD method Vz','Analytical Method')
% % xlabel('Travel time(ms)')
% % ylabel('Amp');
% grid on

%figure 7
Garvin2(2,(150-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVz7(101:end)/10^5-0.5,'r')
temp=conv(diff(src),w)/400;
hold on; plot(temp(101:end)/10^5-0.5,'b')
% % axis([0 3300 -0.5*10^5 0.52*10^5 ])
% legend('Non-Balanced FD method Vz','Analytical Method')
% xlabel('Travel time(ms)')
% ylabel('Amp');
% grid on

%figure 8
Garvin2(2,(200-46)*10/1000)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVz8(101:end)/10^5-0.75,'r')
temp=conv(diff(src),w)/400;
hold on; plot(temp(101:end)/10^5-0.75,'b')
% axis([0 3300 -0.5*10^5 0.52*10^5 ])
% legend('Non-Balanced FD method','Analytical Method')
% xlabel('Travel time(ms)')
% ylabel('Amp');
% grid on

% figure 9
Garvin2(2,2.04)
load('GarvinResult.mat')
hold on;;plot(2.4*4*10^8/400*seismogramVz9(101:end)/10^5-1,'r')
temp=conv(diff(src),w)/400;
hold on; plot(temp(101:end)/10^5-1,'b')
axis([0 3300 -1.5 2.75])
% axis([0 3300 -1.7 3.3])
% legend('Non-Balanced FD method','Analytical Method')
xlabel('Travel time(ms)')
ylabel('Amp');
grid on
