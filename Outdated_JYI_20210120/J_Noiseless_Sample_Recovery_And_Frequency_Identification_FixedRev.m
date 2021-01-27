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
%
% we will s-1 seven frequencies, and make the s-th frequency closer and
% closer to the third one: delta_f=5*10^(-2), 5*10^(-3), 5*10^(-4), 5*10^(-5)

clear all
close all
clc 
warning off

%% Parameters setting up

nc=65;
s=2;
n=2*nc-1;
m=64;
t=0:(n-1);
Omega=randsample(n,m);

a=1;
Q=a*fliplr(diag(ones(nc,1)));

delta_f=1*10^(-3); % in the 1st submission to ISIT 2018 and T-IT, the following are used
% 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001
delta_w=1*10^(-6); % in the 1st submission to ISIT 2018 and T-IT, 10^(-5)
tol=0;

fprintf('\t%8s\t%8s\t%8s\t%8s\n','Nuclear','|Trace|','Tol','Res');

% signal generation
f=zeros(1,s);
f(1)=0.05;
f(2)=f(1)+1/(2*N+1)*10^(-6);
c_mag=[1/2 1/2]';
c_pha=[0 0]';
c=c_mag.*exp(sqrt(-1)*c_pha);
V=exp(1i*2*pi*kron(t',f));
x=V*c; 

% null space criterion
hank_x=hankel(x(1:nc),x(nc:n));
[U_x,D_x,V_x]=svd(hank_x);

U_b=U_x(:,s+1:end);
V_b=V_x(:,s+1:end);
U_x=U_x(:,1:s);
V_x=V_x(:,1:s);

tr=trace(V_x*U_x'*Q);
trabs=abs(tr);
NN=sum(svd(U_b'*Q*V_b));

% signal recovery by EMaC
cvx_solver sdpt3 
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

% tolerance criterion
res_err = norm(u-x)/norm(x);

if res_err<tol && trabs<NN
    fprintf('The signal has been successfully recovered!\n');
    fprintf('The underlying frequencies are being identified...\n');
else 
    fprintf('The signal is not successfully recovered!\n');
end


% post-processing for recovered signal: images for all frequencies in [0,1)
hank_u=hankel(u(1:nc),u(nc:n));
[U_u,D_u,V_u]=svd(hank_u);
U2_u=U_u(:,s+1:end);
w=0:delta_w:1;

J_w=zeros(1,length(w));
for iter_w=1:length(w)

    Atm_w=exp(sqrt(-1)*2*pi*w(iter_w)*[0:nc-1]');
    J_w(iter_w)=(norm(Atm_w,2))/(norm(U2_u'*Atm_w,2));

end

% processing for true signal: images for frequencies in [0,1)
% hank_true=hankel(x(1:nc),x(nc:n));
% [U_true,D_true,V_true]=svd(hank_true);
% U2_true=U_true(:,s+1:end);
% w_true=0:delta_w:1;
% 
% J_w_true=zeros(1,length(w_true));
% for iter_w_true=1:length(w_true)
% 
%     Atm_w_true=exp(sqrt(-1)*2*pi*w_true(iter_w_true)*[0:nc-1]');
%     J_w_true(iter_w_true)=(norm(Atm_w_true,2))/(norm(U2_true'*Atm_w_true,2));
% 
% end

figure 
hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');

figure 
hold on
plot(w,J_w,'-');
xlabel('Frequency f');
ylabel('Imaging function J(f)');
plot(f,zeros(size(f)),'+');
xlim([0.34362 0.34371]);


% figure
% hold on
% stem(w_true,J_w_true);
% xlabel('Frequency f');
% ylabel('Imaging function J(f)');
% stem(f,zeros(size(f)),'+');

% figure
% stem(abs(u-x))
% xlabel('Time instance');
% ylabel('Element difference');

fprintf('The frequency separation is\n');
delta_f
fprintf('The true frequency loation is\n');
f
fprintf('The coefficients are\n');
c
fprintf('The recovery error is\n');
res_err
fprintf('The identified frequency locations are: need look from the figure\n');





