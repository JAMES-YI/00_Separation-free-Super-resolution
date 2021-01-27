% function J_AtomicNormVsHankel()

% This file is to compare performance between atomic norm minimization and Hankel matrix completion 
% for recovering spectrally sparse signal with sparsity R
% 
% The atomic norm minimization will be solved by cvx, and the Hankel matrix
% completion will be first transformed into enhanced matrix completion
% probelm (which is SDP) and then solved by cvx, details can be found in 
% [Chen and Chi 2014, Tang and Bhaskar 2013]
% 
% R: # of exponential terms contained in signal
% Omega: oservation set, here we take it as the whole time domain signal
% value set
% n: # of time instances taken into consideration
% 
% The criterion for success: norm(x-xReconstructed)/norm(x)<10^(-3)
%
% When all signal values are sampled and R takes following values, both methods succeed: nc+1, nc+3, nc+30
% 
% When R=nc+30
%      (1) # of samples n-3: 
% 
% By J. Yi on 10/09/2017, jirong-yi@uiowa.edu


clear all
close all
clc 
warning off

% Parameters setting up

nc=64;
s=nc+3;
n=2*nc-1;
t=0:(n-1);
% Omega=[1:nc-1,nc+1:127];
Omega=randsample(n,80);

% signal generation
f=rand(1,s);
% c=1+10.^(0.5*rand(s,1));
c_mag=1+10.^(0.5*rand(s,1));
c_pha=rand(s,1)*2*pi;
c=c_mag.*exp(sqrt(-1)*c_pha);
V=exp(1i*2*pi*kron(t',f));
x=V*c;

% reconstruction by atomic norm minimization
fprintf('The signal is being reconstructed by atomic norm minimization...\n');
cvx_begin sdp
  variable xAtomic(n) complex;
  variable uAtomic(n) complex;
  variable s(1) complex;
  QAtomic=toeplitz(uAtomic);
  
  minimize real(trace(QAtomic)/n+trace(s))
  subject to 
            [s,xAtomic';xAtomic,QAtomic]>=0;
            xAtomic(Omega)==x(Omega);
cvx_end

resAtomic=norm(x-xAtomic)/norm(x);

% reconstruction by enhanced Hankel matrix completion
fprintf('The signal is being reconstructed by enhanced Hankel matrix completion minimization...\n');
cvx_begin sdp 
    variable Q1(nc,nc) hermitian
    variable Q2(nc,nc) hermitian
    variable xEMaC(n,1) complex 
    QEMaC=[Q1, hankel(xEMaC(1:nc),xEMaC(nc:end)); hankel(xEMaC(1:nc),xEMaC(nc:end))', Q2];
    minimize 0.5*trace(Q1)+0.5*trace(Q2)
    subject to
        QEMaC >= 0,
        xEMaC(Omega) == x(Omega);
cvx_end

resEMaC=norm(x-xEMaC)/norm(x);

% results report
fprintf('The residual in atomic case is %7.2e\n',resAtomic);
fprintf('The residual in Hankel case is %7.2e\n',resEMaC);
% fprintf('True signal is\n');
% x
% fprintf('Signal recovered from atomic norm minimization is\n');
% xAtomic
% fprintf('Signal recovered from enhacned Hankel matrix completion is\n');
% xEMaC


% end