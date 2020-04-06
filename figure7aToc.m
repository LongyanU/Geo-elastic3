
clear;
clc
close all;


load('BalancedTra65Hz.mat')
plotimage(Seismic_Vz(:,46:end-45))
xlabel('x/dx')
ylabel('time(ms)')




load('Nonbalancedf065Hz.mat')
plotimage(Seismic_Vz(:,46:end-45))
xlabel('x/dx')
ylabel('time(ms)')

load('BalancedSchemeM30.mat')
plotimage(Seismic_Vz(:,46:end-45))
xlabel('x/dx')
ylabel('time(ms)')
