% Modified by JYI, 01/26/2021
% - ToDo: documentation

%%

clear all
clc

%% Problem data 
nc=50;
n=2*nc-1;
k=3;
m=30;

%% 1D frequency signal
f0=rand(k,1);
x0=zeros(m,1);
V0=exp(1i*2*pi*kron((0:n-1)',f0'));
c0=rand(k,1)+1i*rand(k,1);
x0=V0*c0;

p=randperm(n);
Omega=p(1:m)';

%% Solve EMaC
tic
cvx_solver sdpt3
cvx_begin sdp quiet
    variable Y(nc,nc) hermitian
    variable Z(nc,nc) hermitian
    variable u(n,1) complex 
    Q1=[Y, hankel(u(1:nc),u(nc:end)); hankel(u(1:nc),u(nc:end))', Z];
    minimize 0.5*trace(Y)+0.5*trace(Z)
    subject to
        Q1 >= 0,
        u(Omega) == x0(Omega);
    
cvx_end
toc
% estimate the frequencies using esprit
res_err = norm(u-x0)/norm(x0)

