% function phase_tran()
% This file is to plot a phase transition of Hankel matriX recovery
% with noiseless measurements with a fixed minimal frequency separation.
%
% - missing one diagonal element
% - s-1 atoms are orthogonal, and the s-th atom has frequency close to that
% of the 1st atom
% 
% Created by JYI, 11/21/2019, based on phase_tran.m
%
clear all; close all; clc; warning off;

%% Parameters setting up
rng(1);
nc = 64; s = 2:1:64;
FName = sprintf('Phase_Tran_Loo_Fig_%s',datestr(now,'yyyymmddHHMM'));

%% Get recovery success rate
% plot phase transition curve
% SuccRateMat = phase_tran_kernel(2:41,88:127,nc,solver); % for phase transition figure
solver = 'HMC'; 
SuccRateMat = phase_tran_kernel_loo(s,nc,solver,100); % for phase transition curve, success rate VS # of atoms
DName = sprintf('%s_Phase_Tran_Loo_Data_%s.mat',...
                 solver,datestr(now,'yyyymmddHHMM'));
save(DName,'SuccRateMat');

figure; hold on;
plot(s,SuccRateMat,'-o');

solver = 'ANM'; 
SuccRateMat = phase_tran_kernel_loo(s,nc,solver,100); % for phase transition curve, success rate VS # of atoms
DName = sprintf('%s_Phase_Tran_Loo_Data_%s.mat',...
                 solver,datestr(now,'yyyymmddHHMM'));
save(DName,'SuccRateMat');

plot(s,SuccRateMat,'-*');
legend('HMC','ANM');
xlabel('# of atoms'); ylabel('Success rate');

%% Plot phase transition
% SuccRateMat = SuccRateMat(svec,mvec);

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3),fig_pos(4)];
print(fig,FName,'-dpdf'); 
close all;


