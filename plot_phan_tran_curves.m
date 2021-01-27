figure; hold on;
load('HMC_Phase_Tran_Loo_Data_201911270309.mat');
plot(2:64,SuccRateMat,'-*');
load('ANM_Phase_Tran_Loo_Data_201911272221.mat')
plot(2:64,SuccRateMat,'-o');
legend('HMC','ANM');
xlabel('# of atoms'); ylabel('Success rate');
xlim([0,70]); ylim([0,1.1])
box on;

% Phase_Tran_Loo_Fig_201911211129.eps

%%
figure; hold on;
load('HMC_phase_tran_orth_loo_data_201911221108.mat');
plot(2:64,succ_prob,'-*');
load('ANM_phase_tran_orth_loo_data_201911291858.mat')
plot(2:64,succ_prob,'-o');
legend('HMC','ANM');
xlabel('# of atoms'); ylabel('Success rate');
xlim([0,70]); ylim([0,1.1])
box on;
% Phase_Tran_Orth_Loo.eps
