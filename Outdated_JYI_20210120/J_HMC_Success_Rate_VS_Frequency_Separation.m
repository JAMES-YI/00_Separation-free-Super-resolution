% This file is to illustrate successful recovery rate for signal with
% s frequency compoents, among which s are randomly generated from [0,1], 
% the smallest distance among frquencies element will be set as Delta_f=1./(2.^[0:1:19]). 
% For each frequency separation, we will randomly do the tests for 20 times, and then plot a
% figure illustrating the recovery performance
% 
% the Hankel matrix has fixed dimension N by N
% 
% the number of measurements is fixed as 2N-4, meaning that we have three
% missing entries
% 
% the number of exponentials is determined by K=w_min/(2*(2N-1-m))-1
% 
% J. Yi on 09/08/2017, jirong-yi@uiowa.edu
% 
% n: # of times instance considered 
% t: time instances considered
% nc: dimension of Hankel matrix
% m: number of measurements
% s: # of exponentials in signal
% delta: smallest frequency separation
% Q: element in null space of sampling operator
% MaxTest: maximal test for each number of measurements
% tol: criterion for determing success of recovery
% ind: index of all time samples not observed
%
%%



function J_HMC_Success_Rate_VS_Frequency_Separation()

clear all
close all
clc 
warning off

%% Parameters setting up
nc=64;

n=2*nc-1;
Full=1:n;
t=0:(n-1);
m=2*nc-4;

while 1
    % s should be fixed
    
    [s,Q,Omega]=Generator(Full,n,m);
    if s>2
        break;
    end
    
end

MaxTest=20;
tol=10^(-5);

Delta_f=[0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001];

fprintf('Separation: \t%8s\t%8s\t%8s\t%8s\t%8s\n','test','Nuclear','|Trace|','Tol','Res');

for iter_f=1:length(Delta_f)
    
    delta_f=Delta_f(iter_f);
    SuccNum=0;
    
    for iter_t=1:MaxTest
        
        % signal generation
        f=randsample([0:8]/9,s-1);
        f(s)=f(max(f)==f)+delta_f;
        c_mag=1+10.^(0.5*rand(s,1));
        c_pha=rand(s,1)*2*pi;
        c=c_mag.*exp(sqrt(-1)*c_pha);
        V=exp(1i*2*pi*kron(t',f));
        x=V*c;
        
        % null space criterion
        hank_x=hankel(x(1:nc),x(nc:n));
        [U_x,~,V_x]=svd(hank_x);

        U_b=U_x(:,s+1:end);
        V_b=V_x(:,s+1:end);
        U_x=U_x(:,1:s);
        V_x=V_x(:,1:s);

        tr=trace(V_x*U_x'*Q);
        trabs=abs(tr);
        NN=sum(svd(U_b'*Q*V_b));
        
        % signal recovery
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

        if res_err<tol && trabs<NN
            SuccNum=SuccNum+1;
        end
        
        % intermediate data
        fprintf('%.2e: \t%d\t%.2e\t%.2e\t%.2e\t%.2e\n',delta_f,iter_t,NN,trabs,tol,res_err);

    end
    
    SuccRate(iter_f)=SuccNum/MaxTest;

end

% graphic illustration
figure
plot(Delta_f,SuccRate,'-*');
xlabel('Frequency separation');
ylabel('Success rate');
axis([0 0.1 0 1.2]);

end

function [s,Q,Omega]=Generator(Full,n,m)

Omega=randsample(n,m);
Unsampled=Full(ismember(Full,Omega)==0);
TopLeft=Unsampled(Unsampled-(n+1)/2<=0);
BottomRight=Unsampled(Unsampled-(n+1)/2>0);
w_min=min(min(TopLeft),n-(max(BottomRight)-1));
s=round(w_min/(2*(n-m))-1);

a=1;
Q=a*fliplr(diag(ones(min(TopLeft),1),(n+1)/2-min(TopLeft)));

end

