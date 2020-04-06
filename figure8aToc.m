
clear;
clc
close all;


load('BalancedTra65Hz.mat')
plotimage(Seismic_Txx(:,46:end-45))
xlabel('x/dx')
ylabel('time(ms)')
title('')



load('Nonbalancedf065Hz.mat')
plotimage(Seismic_Txx(:,46:end-45))
xlabel('x/dx')
ylabel('time(ms)')
title('')


load('BalancedSchemeM30.mat')
plotimage(Seismic_Txx(:,46:end-45))
xlabel('x/dx')
ylabel('time(ms)')
title('')
