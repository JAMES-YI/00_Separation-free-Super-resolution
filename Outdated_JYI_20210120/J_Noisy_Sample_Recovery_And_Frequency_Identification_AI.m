% this file is to illustrate successful recovery for signal with
% s frequency compoents, the Hankel matrix has fixed dimension N by N
%
% 
% J. Yi on 10/10/2017, jirong-yi@uiowa.edu

% n: # of times instance considered 
% t: time instances considered
% nc: dimension of Hankel matrix
% m: # of samples
% s: # of exponentials in signal
% delta: smallest frequency separation
% Q: element in null space of sampling operator
% MaxTest: maximal test for each number of measurements
% tol: criterion for determing success of recovery
% delta_f: minimum frequency separation
% Pn: noise power
%
% Under above setting, the following delta_f can lead to exact
% identification of underlying frequencies: 

clear all
close all
clc 
warning off

%% Parameters setting up

nc=65;
s=8;
n=2*nc-1;
m=64;
t=0:(n-1);
Omega=randsample(n,m);

Pn=0.1;
delta_f=0.001;
tol=10^(-2);

fprintf('\t%8s\t%8s\t%8s\t%8s\n','Nuclear','|Trace|','Tol','Res');

% signal generation
f=zeros(1,s);
f(1:s-1)=rand(1,s-1);
f(s)=f(3)+delta_f;
c_mag=1+10.^(0.5*rand(s,1));
c_pha=rand(s,1)*2*pi;
c=c_mag.*exp(sqrt(-1)*c_pha);
V=exp(1i*2*pi*kron(t',f));
x=V*c; 

% add noise 
nu=randn(n,1)+randn(n,1)*i;
nu=Pn*nu/norm(nu);
xNoisy=x+nu;

% signal recovery by EMaC
cvx_solver sdpt3 
cvx_begin sdp quiet
    variable Y(nc,nc) hermitian
    variable Z(nc,nc) hermitian
    variable u(n,1) complex 
    Q1=[Y, hankel(u(1:nc),u(nc:end)); hankel(u(1:nc),u(nc:end))', Z];
    minimize 0.5*trace(Y)+0.5*trace(Z)
    subject to
        Q1 >= 0,
        u(Omega) == xNoisy(Omega);
cvx_end

% tolerance criterion
res_err = norm(u-x)/norm(x);

% images for all frequencies in [0,1)
hank_u=hankel(u(1:nc),u(nc:n));
[U_u,D_u,V_u]=svd(hank_u);
U2_u=U_u(:,s+1:end);
w=0:delta_f:1;

J_w=zeros(1,length(w));
for iter_w=1:length(w)

    Atm_w=exp(sqrt(-1)*2*pi*w(iter_w)*[0:nc-1]');
    J_w(iter_w)=(norm(Atm_w,2))/(norm(U2_u'*Atm_w,2));

end

% images for underlying frequencies of exponentials
J_f=zeros(1,length(f));
for iter_f=1:length(f)

    Atm_f=exp(sqrt(-1)*2*pi*f(iter_f)*[0:nc-1]');
    J_f(iter_f)=(norm(Atm_f,2))/(norm(U2_u'*Atm_f,2));

end

figure 
hold on
stem(w,J_w);
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');

fprintf('The frequency separation is\n');
delta_f
fprintf('The true frequency loation is\n');
f
fprintf('The coefficients are\n');
c
fprintf('The norm of the noise is\n');
Pn
fprintf('The recovery error is\n');
res_err
fprintf('The identified frequency locations are: need look from the figure\n');









