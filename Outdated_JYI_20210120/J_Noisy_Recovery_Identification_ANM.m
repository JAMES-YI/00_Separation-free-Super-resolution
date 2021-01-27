% this file is to illustrate successful recovery for signal with
% s frequency compoents, the Hankel matrix has fixed dimension N by N
%
% 
% J. Yi on 10/13/2017, jirong-yi@uiowa.edu
% 
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
% we will s-1 seven frequencies, and make the s-th frequency closer and
% closer to the third one

clear all
close all
clc 
warning off

%% Parameters setting up
rng(0);
digits(10);
nc=64;
s=8;
n=2*nc-1;
m=65;
t=0:(n-1);
Omega=randsample(n,m);

Pn=0.1;
delta_f=5*10^(-3); % 0.005 is used in 1st submission to ISIT 2018 and T-IT;
% 5*10^(-6) can be too difficult
delta_w=1*10^(-5);
tol=10^(-2);

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
nu = randn(n,1) + randn(n,1)*i;
nu = Pn*nu/norm(nu);
xNoisy = x+nu;

% signal recovery by EMaC
cvx_begin sdp
  variable u(n,1) complex
  variable ua(n,1) complex
  variable s complex
  Q = toeplitz(ua);
  
  minimize real(trace(Q)/n + trace(s))
  subject to
         [s u';u Q] >= 0;
         u(Omega) == xNoisy(Omega);

cvx_end
clear cvx_*

% tolerance criterion
res_err = norm(u-x)/norm(x);

% images for all frequencies in [0,1)
hank_u=hankel(u(1:nc),u(nc:n));
[U_u,D_u,V_u]=svd(hank_u);
U2_u=U_u(:,s+1:end);

w=0:delta_w:1;
J_w=zeros(1,length(w));
for iter_w=1:length(w)

    Atm_w=exp(sqrt(-1)*2*pi*w(iter_w)*[0:nc-1]');
    J_w(iter_w)=(norm(Atm_w,2))/(norm(U2_u'*Atm_w,2));
    if mod(iter_w,10000) == 0
        fprintf('iter_w: %d/%d\n',iter_w,length(w));
    end

end


% local maxima and frequencies location without dots
figure; 
subplot(2,2,1); hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.00 0.125]);
subplot(2,2,2); hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.125 0.250]);
subplot(2,2,3); hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.250 0.375]);
subplot(2,2,4); hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.375 0.5]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3),fig_pos(4)];
print(fig,'Noisy_ANM_FreIden_P1','-dpdf');
close all;

figure; 
subplot(2,2,1); hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.50 0.625]);
subplot(2,2,2); hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.625 0.750]);
subplot(2,2,3); hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.750 0.875]);
subplot(2,2,4); hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.875 1.00]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3),fig_pos(4)];
print(fig,'Noisy_ANM_FreIden_P2','-dpdf');
close all

figure; hold on;
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.9101 0.9160]);

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3),fig_pos(4)];
print(fig,'Noisy_ANM_FreIden_P3','-dpdf');
close all

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
