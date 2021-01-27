% This file is to illustrate successful recovery rate for signal with
% s frequency compoents which form orthogonal atoms, the Hankel matrix has
% fixed dimension N by N
%
% this file will test 20 times for fixed N to get the successful recovery rate
% 
% J. Yi on 09/06/2017, jirong-yi@uiowa.edu

% n: # of times instance considered 
% t: time instances considered
% nc: dimension of Hankel matrix
% s: # of exponentials in signal
% delta: smallest frequency separation
% Q: element in null space of sampling operator
% MaxTest: maximal test for each number of measurements
% tol: criterion for determing success of recovery
% pctg: vecotr each element of whose is the percentage of measures taken
% 
% Modified by JYI, 01/27/2021
% - all frequencies fall on the grid
% - only the middle sample is missing
% 

%% 

clear all
close all
clc 
warning off

%% Parameters setting up

nc=64;
% s=round(nc-sqrt(log(nc)*nc));
% s=3;
s=20;
n=2*nc-1;
t=0:(n-1);
Omega=[1:nc-1,nc+1:n];

tol=10^(-3);

% fprintf('\t%8s\t%8s\t%8s\t%8s\n','Nuclear','|Trace|','Tol','Res');

% signal generation
f=randsample([0:1/nc:(nc-1)/nc],s);
c_mag=1+10.^(0.5*rand(s,1));
c_pha=rand(s,1)*2*pi;
c=c_mag.*exp(sqrt(-1)*c_pha);
V=exp(1i*2*pi*kron(t',f));
x=V*c; 

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
        u(Omega) == x(Omega);
cvx_end

% tolerance criterion
res_err = norm(u-x)/norm(x);

% visualize frequency locations and imaging function peaks
if res_err<tol
    fprintf('The signal has been successfully recovered!\n');
    fprintf('The underlying frequencies are being identified...\n');
    
        % images for all frequencies in [0,1)
    hank_u=hankel(u(1:nc),u(nc:n));
    [U_u,D_u,V_u]=svd(hank_u);
    U2_u=U_u(:,s+1:end);
    w=0:1/(nc*100):1;
    
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
    plot(w,J_w);
    xlabel('Frequency w');
    ylabel('Imaging J(w)');
    plot(f,zeros(size(f)),'+');
end

% intermediate data
% fprintf('\t%.2e\t%.2e\t%.2e\t%.2e\n',NN,trabs,tol,res_err);




