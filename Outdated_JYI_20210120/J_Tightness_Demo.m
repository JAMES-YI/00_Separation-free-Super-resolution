% This file is to demonstrate the tightness of Theorem 1 for 
% a special case where the sampling set is all but the middle one.
% 
% Created by JYI, 10/24/2019
%
%% 
%% 
clear all
close all
clc 
warning off

%% Parameters setting up
rng(1)
nc=64; 
s=nc/2; % # of atoms
n=2*nc-1; % # of time domain samples in total
t=0:(n-1);
Omega=[1:63, 65:127]';
a=1;
Q=a*fliplr(diag(ones(nc,1)));
delta_w=1*10^(-6); % in the 1st submission to ISIT 2018 and T-IT, 10^(-5)

%% signal generation
f=rand(1,s); %randomly generate frequencies
c_pha = randn(1,s)*2*pi; % coefficient magnitude of atom
c_mag = 1+10.^(0.5*rand(1,s)); % coefficient phase of atom
c=c_mag.*exp(sqrt(-1)*c_pha);
V=exp(1i*2*pi*kron(t',f));
x=V*c'; 

%% signal recovery by EMaC, [Chen et al., 2014]
% cvx_solver sdpt3 
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
%% report results, 
% tolerance criterion
res_err = norm(u-x)/norm(x)

% post-processing for recovered signal: images for all frequencies in [0,1)
% hank_u=hankel(u(1:nc),u(nc:n));
% [U_u,D_u,V_u]=svd(hank_u);
% U2_u=U_u(:,s+1:end);
% w=0:delta_w:1;
% 
% % imaging function for locating frequencies
% J_w=zeros(1,length(w));
% for iter_w=1:length(w)
% 
%     Atm_w=exp(sqrt(-1)*2*pi*w(iter_w)*[0:nc-1]');
%     J_w(iter_w)=(norm(Atm_w,2))/(norm(U2_u'*Atm_w,2));
% 
% end
% 
% % local maxima and frequencies location without dots
% figure; 
% subplot(2,2,1); hold on
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.00 0.125]);
% subplot(2,2,2); hold on
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.125 0.250]);
% subplot(2,2,3); hold on
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.250 0.375]);
% subplot(2,2,4); hold on
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.375 0.5]);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,'Noiseless_Hankel_FreIden_P1','-dpdf');
% close all;
% 
% figure; 
% subplot(2,2,1); hold on
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.50 0.625]);
% subplot(2,2,2); hold on
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.625 0.750]);
% subplot(2,2,3); hold on
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.750 0.875]);
% subplot(2,2,4); hold on
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.875 1.00]);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,'Noiseless_Hankel_FreIden_P2','-dpdf');
% close all;
% 
% figure; hold on;
% plot(w,J_w,'-');
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% plot(f,zeros(size(f)),'+');
% xlim([0.3436 0.3439]);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,'Noiseless_Hankel_FreIden_P3','-dpdf');
% close all;
% 
% fprintf('The frequency separation is\n'); delta_f
% fprintf('The true frequency loation is\n'); f
% fprintf('The coefficients are\n'); c
% fprintf('The recovery error is\n'); res_err
% fprintf('The identified frequency locations are: need look from the figure\n');

