% this file is to illustrate successful recovery rate for signal with
% s frequency compoents which form orthogonal atoms, the Hankel matrix has
% fixed dimension N by N
%
% this file will test 20 times for fixed N to get the successful recovery rate
% 
% J. Yi on 09/04/2017, jirong-yi@uiowa.edu

% n: # of times instance considered 
% t: time instances considered
% nc: dimension of Hankel matrix
% s: # of exponentials in signal
% delta: smallest frequency separation
% Q: element in null space of sampling operator
% MaxTest: maximal test for each number of measurements
% tol: criterion for determing success of recovery
% pctg: vecotr each element of whose is the percentage of measures taken

clear all
close all
clc 
warning off

%% Parameters setting up

nc=64;
s=round(nc-sqrt(log(nc)*nc));
n=2*nc-1;
t=0:(n-1);
Omega=[1:63,65:n];

a=1;
Q=a*fliplr(diag(ones(nc,1)));

MaxTest=20;
tol=10^(-3);

fprintf('Test: \t%8s\t%8s\t%8s\t%8s\n','Nuclear','|Trace|','Tol','Res');

NumSucc=0;
for iter_m=1:MaxTest
        
    % signal generation
    f=randsample([0:1/nc:(nc-1)/nc],s);
    c_mag=1+10.^(0.5*rand(s,1));
    c_pha=rand(s,1)*2*pi;
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
        NumSucc=NumSucc+1;
    end

    % intermediate data
    fprintf('%4d: \t%.2e\t%.2e\t%.2e\t%.2e\n',iter_m,NN,trabs,tol,res_err);
    
end

% successful recovery rate
RateSucc=NumSucc/20


