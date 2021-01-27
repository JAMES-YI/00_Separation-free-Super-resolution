% This file is to demonstrate the tightness for successful recovery when
% the frequency atoms are orthogonal.
% 
% Created by JYI, 11/11/2019
%
%
%%

clear all; close all; clc; warning off;

%% Parameters setting up
rng(1);

nc = 64; n = 2*nc-1; fVec = 0:1/nc:1;
t = 0:(n-1); nTrial = 20;
succ_prob = [];
m = 95;
DName = sprintf('phase_tran_orth_data_m%d_%s.mat',...
    m,datestr(now,'yyyymmddHHMM'));

for s=2:1:64
    
    succ_count = 0;
    
    for i=1:nTrial
         
        f = randsample(fVec,s);
        Omega = randsample(n,m);
        %% Signal generation

        c_mag = 1+10.^(0.5*rand(s,1));
        c_pha = 2*pi*rand(s,1);
        c = c_mag.*exp(sqrt(-1)*c_pha);
        V = exp(1i*2*pi*kron(t',f));
        x = V*c; 


        %% Signal recovery by EMaC (Hankel matrix recovery approach), [Chen et al., 2014]

        cvx_begin sdp quiet
            variable Y(nc,nc) hermitian
            variable Z(nc,nc) hermitian
            variable u(n,1) complex 
            Q1=[Y, hankel(u(1:nc),u(nc:end)); hankel(u(1:nc),u(nc:end))', Z];
            minimize 0.5*trace(Y)+0.5*trace(Z)
            subject to
                Q1 >= 0,
                u(Omega) == x(Omega);
        cvx_end
        clear cvx_*  
        
        res_err = norm(u-x)/norm(x);
        fprintf('s: %d, trial: %d/%d, residual: %4.2e\n',s,i,nTrial,res_err);
        
        if res_err<10^(-6)
            succ_count = succ_count+1;
        end
    end
    
    succ_prob = [succ_prob,succ_count/nTrial]
end
save(DName,'succ_prob');



%% Report results
data95 = load('phase_tran_orth_data_m95_201911181409.mat');
data65 = load('phase_tran_orth_data_m65_201911171033.mat');
data35 = load('phase_tran_orth_data_m35_201911151418.mat');

figure; hold on;
plot(2:1:64,data95.succ_prob,'-*');
plot(2:1:64,data65.succ_prob,'-o');
plot(2:1:64,data35.succ_prob,'-<');

xlabel('# of atoms'); ylabel('Success rate');
legend('Measurements = 95','Measurements = 65','Measurements = 35');

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3),fig_pos(4)];
print(fig,'Tightness_Orth','-dpdf');
close all;

