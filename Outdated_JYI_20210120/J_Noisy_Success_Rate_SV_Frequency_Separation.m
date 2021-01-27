% this file is to illustrate successful recovery rate for signal with
% s frequency compoents randomly generated from [0,1), the Hankel matrix has
% fixed dimension N by N
%
% this file will test 100 times for fixed N to get the successful recovery rate
% 
% J. Yi on 10/10/2017, jirong-yi@uiowa.edu

% n: # of times instance considered 
% t: time instances considered
% nc: dimension of Hankel matrix
% m: # of measurements
% nu: noise with ||nu||_2
% Pn: noise power level
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
s=8;
m=64;
n=2*nc-1;
t=0:(n-1);
Pn=0.001;
Delta_f=[0.05,0.03,0.01,0.005,0.003,0.001,...
         0.0005,0.0003,0.0001,0.00005,0.00003,0.00001];

MaxTest=100;
tol=10^(-2);

filename='J_Noisy_Success_Rate_VS_Frequency_Separation';
col_name={'Sepa','TestNumb','Tol','Residual','Succ'};
rowInd=1;
Range=sprintf('A%d',rowInd);
xlswrite(filename,col_name,'sheet',Range);

for iter_f=1:length(Delta_f)

    delta_f=Delta_f(iter_f);
    fprintf('Frequency separation is %.2e\n',delta_f);
    fprintf('Test: \t%8s\t%8s\t%8s\n','Tol','Res','Success');
    NumSucc=0;

    for iter_m=1:MaxTest
        
        % signal generation
        f=zeros(1,s);
        f(1:s-1)=randsample([0:1/10:(10-1)/10],s-1);
        f(s)=f(s-4)+delta_f;
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
        Omega=randsample(1:1:n,m);
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

        if res_err<tol 
            NumSucc=NumSucc+1;
            Succ=1;
        else
            Succ=0;
        end

        % result report
        fprintf('%4d: \t%.2e\t%.2e\t%.2e\n',iter_m,tol,res_err,Succ);
        
        rowInd=rowInd+1;
        Entry=[delta_f,iter_m,tol,res_err,Succ];
        Range=sprintf('A%d',rowInd);
        xlswrite(filename,Entry,'sheet',Range);
        
    end
    
    RateSucc(iter_f)=NumSucc/100;

end

plot(Delta_f,RateSucc,'-*');
xlabel('Frequency separation');
ylabel('Succeeful recovery rate');



