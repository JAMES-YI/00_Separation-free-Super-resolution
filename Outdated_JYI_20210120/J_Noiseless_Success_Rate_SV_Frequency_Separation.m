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
Delta_f=[0.05,0.03,0.01,0.005,0.003,0.001,...
         0.0005,0.0003,0.0001,0.00005,0.00003,0.00001];

a=1;
Q=a*fliplr(diag(ones(nc,1)));

MaxTest=100;
tol=10^(-5);

filename='J_Noiseless_Success_Rate_VS_Frequency_Separation';
col_name={'Sepa','TestNumb','NuNorm','AbsTrace','Tol','Residual','Succ'};
rowInd=1;
Range=sprintf('A%d',rowInd);
xlswrite(filename,col_name,'sheet',Range);

for iter_f=1:length(Delta_f)

    delta_f=Delta_f(iter_f);
    fprintf('Frequency separation is %.2e\n',delta_f);
    fprintf('Test: \t%8s\t%8s\t%8s\t%8s\t%8s\n','Nuclear','|Trace|','Tol','Res','Success');
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
                u(Omega) == x(Omega);
        cvx_end

        % tolerance criterion
        res_err = norm(u-x)/norm(x);

        if res_err<tol && trabs<NN
            NumSucc=NumSucc+1;
            Succ=1;
        else
            Succ=0;
        end

        % result report
        fprintf('%4d: \t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n',iter_m,NN,trabs,tol,res_err,Succ);
        
        rowInd=rowInd+1;
        Entry=[delta_f,iter_m,NN,trabs,tol,res_err,Succ];
        Range=sprintf('A%d',rowInd);
        xlswrite(filename,Entry,'sheet',Range);
        
    end
    
    RateSucc(iter_f)=NumSucc/100;

end

plot(Delta_f,RateSucc,'-*');
xlabel('Frequency separation');
ylabel('Succeeful recovery rate');



