% This file is for running the experiments on argon cluster using high
% performance computation in a parallel way.
%
% Created by JYI, 10/28/2019
%
%%
clear all; close all; clc; warning off;

%% Parameters setting up

nc = 64; 
solver = 'ANM'; 

parpool('open',20);
parfor i=1:20
    
    SuccRateMat = phase_tran_kernel(2*i:2i+1,88:127,nc,solver,20); % for phase transition curve, success rate VS # of atoms

end
parpool('close');

DName = sprintf('%s_Phase_Tran_Data_%s.mat',...
                 solver,datestr(now,'yyyymmddHHMM'));
save(DName,'SuccRateMat');
