% This file is to demonstrate the performance of Hankel matrix recovery
% approach in two situations:
% (1) with noiseless undersampled observations.
% (2) with noisy undersampled observations.
% 
% - nc: dimension of Hankel matrix
% - n: # of times instance considered 
% - m: # of observed samples
% - Omega: index set of observed samples
% - delta_f: minimum frequency separation
% - delta_w: precision when using MUSIC for identifying frequencies
% - s: # of exponentials in signal
% - t: time indices considered
% - f: frequencies
% - V: atoms
% - c: coefficients
% - x: ground truth signal
% - u: recovered signal
%
% - (s-1) out of s frequencies are prefixed, and the rest one is taken to
%   be close to the 3rd frequency
% - the CVX solver is needed
% - change delta_f to control the minimal frequency separation
%
% Created by JYI, 10/24/2019, jirong-yi@uiowa.edu, based on the version created on 10/10/2017 by JYI.

clear all; close all; clc; warning off;

%% Parameters setting up
rng(1)

nc = 64; n = 2*nc-1; 
m=65; Omega=randsample(n,m); Jcomp = setdiff(1:n, Omega);

% delta_f = 5*10^(-3); % in the 1st submission to ISIT 2018 and T-IT, the following are used
% delta_f_Vec = [0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001];
delta_w = 1*10^(-6); % in the 1st submission to ISIT 2018 and T-IT, 10^(-5)
s = 8; t = 0:(n-1);

% Observation type
obstype = 'Noisy'; 
if strcmp(obstype,'Noisy')
    Pn = 1.5;
    delta_f_Vec = 5*10^(-3);
elseif strcmp(obstype,'Noiseless')
    JYI = 0;
    delta_f_Vec = [0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001];
end

% Choose solver
solver_Vec = {'ANM','HMC'};
% solver = 'ANM'; % or HMC, ANM
% figname = sprintf('%s_%s_sep_5_e-3_log',solver,obstype);
%% Signal generation
for i=1:2
    solver = solver_Vec(i);
    
    for delta_f=delta_f_Vec
        
        
        sep = num2str(delta_f);
        figname = sprintf('%s_%s_sep_%s_log',...
            solver{1},obstype,sep(3:end));
        
        f = zeros(1,s);
        f(1:s-1) = [0.3923, 0.9988, 0.3437, 0.9086, 0.6977, 0.0298, 0.4813];
        f(s)=f(3)+delta_f;
        c_mag = [3.1800, 2.5894, 2.1941, 2.9080, 3.9831, 4.0175, 4.1259, 3.6182]';
        c_pha = [4.1097, 5.4612, 5.4272, 4.7873, 1.0384, 0.4994, 3.1975, 0.5846]';
        c = c_mag.*exp(sqrt(-1)*c_pha);
        V = exp(1i*2*pi*kron(t',f));
        x = V*c; 

        if strcmp(obstype,'Noisy')

            nu = randn(n,1) + randn(n,1)*sqrt(-1);
            nu = Pn*nu/norm(nu);
            x = x+nu;
        end

        %% Signal recovery by EMaC (Hankel matrix recovery approach), [Chen et al., 2014]
        % or ANM, [Tang et al., 2013]

        if strcmp(solver,'HMC')
            cvx_begin sdp 
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

            % Identify frequency via MUSIC
            hank_u = hankel(u(1:nc),u(nc:n));
            [U_u,D_u,V_u] = svd(hank_u);
            U2_u = U_u(:,s+1:end);
            w = 0:delta_w:1;

            Atm = exp(sqrt(-1)*2*pi*kron((0:nc-1)',w));
            Denom = vecnorm(U2_u'*Atm,2,1);
            Nume = vecnorm(Atm,2,1);
            J_w = log10(Nume ./ Denom);

            % Visualize frequency identification performance
            figure; 
            subplot(2,1,1); hold on;
            plot(w,J_w,'-'); plot(f,zeros(size(f)),'+');
            xlabel('Frequency'); ylabel('Logarithm of imaging function');
            box on;
            subplot(2,1,2); hold on;
            plot(w,J_w,'-'); plot(f,zeros(size(f)),'+'); xlim([0.3426 0.3740]);
            xlabel('Frequency'); ylabel('Logarithm of imaging function');
            box on;


            fig = gcf;
            fig.PaperPositionMode = 'auto';
            fig_pos = fig.PaperPosition;
            fig.PaperSize = [fig_pos(3),fig_pos(4)];
            print(fig,figname,'-dpdf');
            close all;
        elseif strcmp(solver,'ANM')

            % solve primal problem of atomic norm minimization
            cvx_begin sdp quiet
            variable u(n,1) complex
            variable ua(n,1) complex
            variable sv complex
            Q = toeplitz(ua);

            minimize real(trace(Q)/n + trace(s))
            subject to
                 [sv u';u Q] >= 0;
                 u(Omega) == x(Omega);

            cvx_end
            clear cvx_*

            % solve dual problem of ANM
            cvx_begin quiet
                variable q(n) complex;
                variable H(n,n) hermitian;

                maximize (real(x'*q))

                subject to
                    [H q; q' 1] == semidefinite(n+1,n+1);
                    trace(H) == 1;
                    for j = 1:n-1
                        sum(diag(H,j)) == 0;
                    end
                    q(Jcomp) == 0;
            cvx_end

            % Identify frequency via dual polynomial
            N = 2^16;
            freqgrid = (0:N-1)'/N;
            af = exp(1j*2*pi*t'*freqgrid');
            Qf = af'*q;

            figure; 
            subplot(2,1,1); hold on;
            plot(f,zeros(1,s),'+'); plot(freqgrid,abs(Qf),'-');
            xlabel('Frequency'); ylabel('Magnitude of dual polynomial');
            ylim([0,1.1]); box on;
            
            subplot(2,1,2); hold on;
            plot(f,zeros(1,s),'+'); plot(freqgrid,abs(Qf),'-'); xlim([0.3426 0.3740]);
            xlabel('Frequency'); ylabel('Magnitude of dual polynomial');
            ylim([0,1.1]); box on;

            fig = gcf;
            fig.PaperPositionMode = 'auto';
            fig_pos = fig.PaperPosition;
            fig.PaperSize = [fig_pos(3),fig_pos(4)];
            print(fig,figname,'-dpdf');
            close all;
        end


        %% Report results, 

        res_err = norm(u-x)/norm(x);
        err = norm(u-x); 
        fprintf('Relative reconstruction error by %s: %.4e, reconstruction error: %.4e\n',solver{1},res_err,err)

    end
end
% f = zeros(1,s);
% f(1:s-1) = [0.3923, 0.9988, 0.3437, 0.9086, 0.6977, 0.0298, 0.4813];
% f(s)=f(3)+delta_f;
% c_mag = [3.1800, 2.5894, 2.1941, 2.9080, 3.9831, 4.0175, 4.1259, 3.6182]';
% c_pha = [4.1097, 5.4612, 5.4272, 4.7873, 1.0384, 0.4994, 3.1975, 0.5846]';
% c = c_mag.*exp(sqrt(-1)*c_pha);
% V = exp(1i*2*pi*kron(t',f));
% x = V*c; 
% 
% if strcmp(obstype,'Noisy')
%     
%     nu = randn(n,1) + randn(n,1)*sqrt(-1);
%     nu = Pn*nu/norm(nu);
%     x = x+nu;
% end
% 
% %% Signal recovery by EMaC (Hankel matrix recovery approach), [Chen et al., 2014]
% % or ANM, [Tang et al., 2013]
% 
% if strcmp(solver,'HMC')
%     cvx_begin sdp 
%         variable Y(nc,nc) hermitian
%         variable Z(nc,nc) hermitian
%         variable u(n,1) complex 
%         Q1=[Y, hankel(u(1:nc),u(nc:end)); hankel(u(1:nc),u(nc:end))', Z];
%         minimize 0.5*trace(Y)+0.5*trace(Z)
%         subject to
%             Q1 >= 0,
%             u(Omega) == x(Omega);
%     cvx_end
%     clear cvx_*
%     
%     % Identify frequency via MUSIC
%     hank_u = hankel(u(1:nc),u(nc:n));
%     [U_u,D_u,V_u] = svd(hank_u);
%     U2_u = U_u(:,s+1:end);
%     w = 0:delta_w:1;
% 
%     Atm = exp(sqrt(-1)*2*pi*kron((0:nc-1)',w));
%     Denom = vecnorm(U2_u'*Atm,2,1);
%     Nume = vecnorm(Atm,2,1);
%     J_w = log10(Nume ./ Denom);
%     
%     % Visualize frequency identification performance
%     figure; 
%     subplot(2,1,1); hold on;
%     plot(w,J_w,'-'); plot(f,zeros(size(f)),'+');
%     xlabel('Frequency'); ylabel('Logarithm of imaging function');
%     subplot(2,1,2); hold on;
%     plot(w,J_w,'-'); plot(f,zeros(size(f)),'+'); xlim([0.3426 0.3740]);
%     xlabel('Frequency'); ylabel('Logarithm of imaging function');
%     
% 
%     fig = gcf;
%     fig.PaperPositionMode = 'auto';
%     fig_pos = fig.PaperPosition;
%     fig.PaperSize = [fig_pos(3),fig_pos(4)];
%     print(fig,figname,'-dpdf');
%     close all;
% elseif strcmp(solver,'ANM')
%     
%     % solve primal problem of atomic norm minimization
%     cvx_begin sdp
%     variable u(n,1) complex
%     variable ua(n,1) complex
%     variable sv complex
%     Q = toeplitz(ua);
% 
%     minimize real(trace(Q)/n + trace(s))
%     subject to
%          [sv u';u Q] >= 0;
%          u(Omega) == x(Omega);
% 
%     cvx_end
%     clear cvx_*
%     
%     % solve dual problem of ANM
%     cvx_begin
%         variable q(n) complex;
%         variable H(n,n) hermitian;
% 
%         maximize (real(x'*q))
% 
%         subject to
%             [H q; q' 1] == semidefinite(n+1,n+1);
%             trace(H) == 1;
%             for j = 1:n-1
%                 sum(diag(H,j)) == 0;
%             end
%             q(Jcomp) == 0;
%     cvx_end
%     
%     % Identify frequency via dual polynomial
%     N = 2^16;
%     freqgrid = (0:N-1)'/N;
%     af = exp(1j*2*pi*t'*freqgrid');
%     Qf = af'*q;
%     
%     figure; 
%     subplot(2,1,1); hold on;
%     plot(f,zeros(1,s),'+'); plot(freqgrid,abs(Qf),'-');
%     xlabel('Frequency'); ylabel('Magnitude of dual polynomial');
%     subplot(2,1,2); hold on;
%     plot(f,zeros(1,s),'+'); plot(freqgrid,abs(Qf),'-'); xlim([0.3426 0.3740]);
%     xlabel('Frequency'); ylabel('Magnitude of dual polynomial');
% 
%     fig = gcf;
%     fig.PaperPositionMode = 'auto';
%     fig_pos = fig.PaperPosition;
%     fig.PaperSize = [fig_pos(3),fig_pos(4)];
%     print(fig,figname,'-dpdf');
%     close all;
% end
%         
% 
% %% Report results, 
% 
% res_err = norm(u-x)/norm(x);
% err = norm(u-x); 
% fprintf('Relative reconstruction error by %s: %.4e, reconstruction error: %.4e\n',solver,res_err,err)
%% Visual results: set 0


%% Visual results: set 1
% ZoomName = sprintf('%s_%s_Sep_1_e-4',solver,obstype);
% 
% Split = 1; clear Names; 
% Names.F1 = sprintf('%s_%s_Sep_1_e-4_P1',solver,obstype); 
% Names.F2 = sprintf('%s_%s_Sep_1_e-4_P2',solver,obstype);
%Split = 0; clear Names; Names = sprintf('%s_%s_Sep_1_e-4_F',solver,obstype);

% fre_loc_vis(Split,Names,w,J_w,f)

% figure; hold on;
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.3426 0.3440]);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,ZoomName,'-dpdf');
% close all;

%% Phase transition
% hmc = load('HMC_Phase_Tran_Data_201910261221.mat');
% anm = load('ANM_Phase_Tran_Data_201910260032.mat');
% figure; hold on;
% plot(hmc.SuccRateMat(:,65),'-o');
% plot(anm.SuccRateMat(:,65),'-*');
% legend('HMC','ANM');
% xlabel('# of atoms');
% ylabel('Success rate');
% 
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,'Phase_Tran_Fig','-dpdf');




