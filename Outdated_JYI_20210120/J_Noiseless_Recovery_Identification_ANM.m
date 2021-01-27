%% This file is to illustrate successful recovery for signal with
% s frequency compoents, the Hankel matrix has fixed dimension N by N
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
%
% we will s-1 seven frequencies, and make the s-th frequency closer and
% closer to the third one: delta_f
%
% the CVX solver is needed.
%
% J. Yi on 10/10/2017, jirong-yi@uiowa.edu

%% 
clear all
close all
clc 
warning off

%% Parameters setting up
rng(1)
nc=64;
s=8;
n=2*nc-1;
m=65;
t=0:(n-1);
Omega=randsample(n,m);
a=1;
Q=a*fliplr(diag(ones(nc,1)));
delta_f=10^(-4); 
delta_w=1*10^(-6); 
tol=0;

%% signal generation
f=zeros(1,s);
f(1:s-1)=[0.3923, 0.9988, 0.3437, 0.9086, 0.6977, 0.0298, 0.4813];
f(s)=f(3)+delta_f;
c_mag=[3.1800, 2.5894, 2.1941, 2.9080, 3.9831, 4.0175, 4.1259, 3.6182]';
c_pha=[4.1097, 5.4612, 5.4272, 4.7873, 1.0384, 0.4994, 3.1975, 0.5846]';
c=c_mag.*exp(sqrt(-1)*c_pha);
V=exp(1i*2*pi*kron(t',f));
x=V*c; 

%% signal recovery by ANM, [Tang et al., 2013]
% solution is u

cvx_begin sdp
  variable u(n,1) complex
  variable ua(n,1) complex
  variable s complex
  Q = toeplitz(ua);
  
  minimize real(trace(Q)/n + trace(s))
  subject to
         [s u';u Q] >= 0;
         u(Omega) == x(Omega);

cvx_end
clear cvx_*
%% report results, 
% tolerance criterion
res_err = norm(u-x)/norm(x);

% post-processing for recovered signal: images for all frequencies in [0,1)
hank_u=hankel(u(1:nc),u(nc:n));
[U_u,D_u,V_u]=svd(hank_u);
U2_u=U_u(:,s+1:end);
w=0:delta_w:1;

% imaging function for locating frequencies
J_w=zeros(1,length(w));
for iter_w=1:length(w)

    Atm_w=exp(sqrt(-1)*2*pi*w(iter_w)*[0:nc-1]');
    J_w(iter_w)=(norm(Atm_w,2))/(norm(U2_u'*Atm_w,2));

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
print(fig,'Noiseless_ANM_FreIden_P1','-dpdf');
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
print(fig,'Noiseless_ANM_FreIden_P2','-dpdf');
close all;

figure; hold on;
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.3436 0.3439]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3),fig_pos(4)];
print(fig,'Noiseless_ANM_FreIden_P3','-dpdf');
close all;

fprintf('The frequency separation is\n'); delta_f
fprintf('The true frequency loation is\n'); f
fprintf('The coefficients are\n'); c
fprintf('The recovery error is\n'); res_err
fprintf('The identified frequency locations are: need look from the figure\n');





