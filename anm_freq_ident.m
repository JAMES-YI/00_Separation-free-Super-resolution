% This file is to demonstrate the peformance of recovering spectrally sparse signal 
% and identifying underlying frequencies via atomic norm minimization in
% dual form & Hankel matrix completion with MUSIC.
%
% References
% [1] [Tang et al., 2013]
% [2] [Xu et al., 2017]
%
% Codes adapted from Kumar Vijay Mishra and Myung Cho.
% 
% Created by JYI, 11/19/2019
%
% Modified by JYI, 01/24/2021
% - ToDo: documentation

%%
clear all; close all; clc; warning off;

%% Parameters setting up
rng(1);

nc = 64; n = 2*nc-1; t = 0:(n-1); 
s = 16; 
m = 95; 
del_f = 10^(-5); % minimal frequency separation 
del_w = 10^(-6); 
zoomxrange = [0.2791,0.2792];
anmfigname = sprintf('ANM_Freq_Iden_Sep_e-5');
hmcfigname = sprintf('HMC_Freq_Iden_Sep_e-5');

Omega = randsample(n,m); Jcomp = setdiff(1:n, Omega);
%% Signal generation
f = zeros(1,s); % frequencies
f(1:s-1) = rand(1,s-1);
f(s) = f(3) + del_f;
if f(s)>1
    f(s) = f(3) - destel_f;
end

% fprintf('%.8f,%.8f\n',f(3),f(end));

c_mag = 1+10.^(0.5*rand(s,1));
c_pha = 2*pi*rand(s,1);
c = c_mag.*exp(sqrt(-1)*c_pha); % coefficient of each atom
V = exp(1i*2*pi*kron(t',f)); % each column is an atom
x = V*c; % original signal

%% Signal recovery by ANM [Tang et al., 2013]

% solve primal problem of ANM, eq(II.6) in [Tang et al., 2013]
%     minimize_x,u,t (1/(2*length(x)) * trace(Toep(u)) + (1/2) * t
%     subject to     [Toep(u), x; x', t] >=0
%                    x_j = x_j^*, j in T
%     - x^*_j: observation of the j-th element of the ground truth
%     - x: reconstructed signal
%     - u: a complex vector
%     - Toep(u): Toeplitz matrix generated from u, with first column being u, and first row being u^*;
%       the function toeplitz(u) in MATLAB will generate a Toeplitz matrix
%       with first column being u^*, and first row being u; the (u)_1 and
%       (u^*)_1 may not be the same, and we take (u)_1 as (toeplitz(u))_11;
%       this also results in a different form for implementation as [s u';u Q] >= 0;


% cvx_begin sdp quiet
%     variable u(n,1) complex
%     variable ua(n,1) complex
%     variable vs complex
%     Q = toeplitz(ua);
% 
%     minimize real(trace(Q)/n + trace(s))
%     subject to
%          [vs u';u Q] >= 0;
%          u(Omega) == x(Omega);
% 
% cvx_end
% clear cvx_*
    
% solve dual problem of ANM, eq(II.13) in [Tang et al., 2013]
%     maximize_q,H <q_Omega,x_Omega^*>_R
%     subject to    [H,-q;-q^*,1] >= 0
%                   sum_k=1^{|Omega|-j} H_k,k+j = 1, j=0
%                   sum_k=1^{|Omega|-j} H_k,k+j = 0, j=1,2,...,|Omega|-1
%                   H is Hermitian
%                   q_[n]\Omega = 0

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

%% Identify frequencies from dual polynomial
N = 2^16;
freqgrid = (0:N-1)'/N;
af = exp(1j*2*pi*t'*freqgrid');
Qf = af'*q;

figure; 
subplot(2,1,1); hold on;
plot(f,zeros(1,s),'+');
plot(freqgrid,abs(Qf),'-');
xlabel('Frequency'); ylabel('Magnitude of dual polynomial');
subplot(2,1,2); hold on;
plot(f,zeros(1,s),'+'); plot(freqgrid,abs(Qf),'-'); xlim(zoomxrange)
xlabel('Frequency'); ylabel('Magnitude of dual polynomial');

% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,anmfigname,'-dpdf');
% close all;

%% Solve Hankel matrix recovery, [Xu et al., 2017]
cvx_begin sdp quiet
    variable Y(nc,nc) hermitian
    variable Z(nc,nc) hermitian
    variable uhmc(n,1) complex 
    Q1=[Y, hankel(uhmc(1:nc),uhmc(nc:end)); hankel(uhmc(1:nc),uhmc(nc:end))', Z];
    minimize 0.5*trace(Y)+0.5*trace(Z)
    subject to
        Q1 >= 0,
        uhmc(Omega) == x(Omega);
cvx_end
clear cvx_*

%% Identify frequencies via MUSIC
hank_u = hankel(uhmc(1:nc),u(nc:n));
[U_u,D_u,V_u] = svd(hank_u);
U2_u = U_u(:,s+1:end);
w = 0:del_w:1;

Atm = exp(sqrt(-1)*2*pi*kron((0:nc-1)',w));
Denom = vecnorm(U2_u'*Atm,2,1);
Nume = vecnorm(Atm,2,1);
J_w = log10(Nume ./ Denom);

figure; 
subplot(2,1,1); hold on;
plot(w,J_w,'-'); plot(f,zeros(size(f)),'+');
xlabel('Frequency'); ylabel('Logarithm of imaging function');
subplot(2,1,2); hold on;
plot(w,J_w,'-'); plot(f,zeros(size(f)),'+'); xlim(zoomxrange);
xlabel('Frequency'); ylabel('Logarithm of imaging function');

% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,hmcfigname,'-dpdf');
% close all;

%% Quantitative results
resanm = norm(u-x) / norm(x);
reshmc = norm(uhmc-x) / norm(x);
fprintf('Relative error, ANM: %8.4e, HMC: %8.4e\n',resanm,reshmc);




