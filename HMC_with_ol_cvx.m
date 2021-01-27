clear all
clc

%% Problem data 
nc=50;
n=2*nc-1;
k=3;
m=70;

%% 1D frequency signal
f0=rand(k,1);
x0=zeros(m,1);
V0=exp(1i*2*pi*kron((0:n-1)',f0'));
c0=rand(k,1)+1i*rand(k,1);
x0=V0*c0;

p=randperm(n);
Omega=p(1:m)';

% outlier
no=3;
ol=100*(rand(no,1)+1i*rand(no,1));
po=randperm(m);
o0=zeros(n,1);
o0(Omega(po(1:no)'))=ol; % outlier is added to observations over the only omega set.

% final signal
y0=x0+o0;

%% Solve EMaC with sparse outlier
lambda=2;
tic
cvx_solver sdpt3
cvx_begin sdp quiet
    variable Y2(nc,nc) hermitian
    variable Z2(nc,nc) hermitian
    variable u2(n,1) complex 
    variable o2(n,1) complex
    Q2=[Y2, hankel(u2(1:nc),u2(nc:end)); hankel(u2(1:nc),u2(nc:end))', Z2];
    y2=u2+o2;
    minimize 0.5*trace(Y2)+0.5*trace(Z2)+lambda*norm(o2,1)
    subject to
        Q2 >= 0,
        y2(Omega) == y0(Omega);
cvx_end
toc
% estimate the frequencies using esprit
res_err2 = norm(y2-y0)/norm(y0)