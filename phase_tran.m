% function phase_tran()
% This file is to plot a phase transition of Hankel matriX recovery
% with noiseless measurements with a fixed minimal frequency separation.
%
% - % (2:64,65) is used for obtaining
%      HMC_Phase_Tran_Data_201910261221.mat, 
%      ANM_Phase_Tran_Data_201910260032.mat
%
% Created by JYI, 10/25/2019.
%
clear all; close all; clc; warning off;

%% Parameters setting up
rng(1);
nc = 64; 
s = 2:1:64; m= 65; % 35, 65, 95, 125
FName = sprintf('Phase_Tran_Fig_m%d_%s',m,datestr(now,'yyyymmddHHMM'));


%% Get recovery success rate
% plot phase transition curve
% SuccRateMat = phase_tran_kernel(2:41,88:127,nc,solver); % for phase transition figure
solver = 'HMC'; 
SuccRateMat = phase_tran_kernel(s,m,nc,solver,100); % for phase transition curve, success rate VS # of atoms

DName = sprintf('%s_Phase_Tran_Data_m%d_%s.mat',...
                 solver,m,datestr(now,'yyyymmddHHMM'));
save(DName,'SuccRateMat');

figure; hold on;
plot(s,SuccRateMat(:,1),'-o');

solver = 'ANM'; 
SuccRateMat = phase_tran_kernel(s,m,nc,solver,20); % for phase transition curve, success rate VS # of atoms
DName = sprintf('%s_Phase_Tran_Data_m%d_%s.mat',...
                 solver,m,datestr(now,'yyyymmddHHMM'));
save(DName,'SuccRateMat');

plot(s,SuccRateMat(:,1),'-*');
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


