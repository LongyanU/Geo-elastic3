
clear;
clc
close all;


load('BalancedTra65Hz.mat')
% plotimage(Seismic_Vx)
a=Seismic_Vz;


% figure;plot(Seismic_Vx(:,267))
load('Nonbalancedf065Hz.mat')
b=Seismic_Vz;

 load('BalancedSchemeM30.mat')
c=Seismic_Vz;

figure;
% plot(a(:,267),'b','linewidth',1)
% hold on;plot(b(:,267),'k','linewidth',1)
% hold on;plot(c(:,267),'r','linewidth',1)
hold on;plot(a(:,267)-c(:,267),'b','linewidth',1)
hold on;plot(b(:,267)-c(:,267),'k','linewidth',1)

legend('Difference between balanced SGFD and reference','Difference between non-balanced SGFD and reference')
grid on
xlabel('travel time(ms)')
box on
ylabel('Difference of Vz(m/s)')
% axis([0 4000 -6*10^-11  6*10^-11])
axis([0 4000 -6.6*10^-10  9.2*10^-10])
% figure;plot(a(:,300))
% hold on;plot(b(:,300),'k')
% hold on;plot(c(:,300),'r')
% hold on;plot(c(:,300)-b(:,300)-1*10^-9,'m')
% hold on;plot(a(:,300)-b(:,300)-2*10^-9,'m')
% legend('Balanced FD Scheme M=7','Balanced FD Scheme M=22','Non-Balanced FD Scheme M=7')
% 
% figure;plot(a(:,500))
% hold on;plot(b(:,500),'k')
% hold on;plot(c(:,500),'r')
% hold on;plot(c(:,500)-b(:,500)-1*10^-9,'m')
% hold on;plot(a(:,500)-b(:,500)-2*10^-9,'m')
% legend('Balanced FD Scheme M=7','Balanced FD Scheme M=22','Non-Balanced FD Scheme M=7')