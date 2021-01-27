% this file is to illustrate successful recovery rate for signal with
% s frequency compoents, and s-1 of them form orthogonal atoms, while the
% rest one frequency is very close to one of the s-1 frequencies with
% separation delta
%
% specifically, this file will give performance illustration of the
% relation between successful recovery rate and number of measurements. the separation delta is fixed
% 
% J. Yi on 08/04/2017, jirong-yi@uiowa.edu

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
n=127;
t=0:(n-1);
nc=(n+1)/2;
s=round(0.05*n);
a=1;
delta=1/512;
Q=a*fliplr(diag(1,-(nc-1)));

MaxTest=10;
tol=10^(-6);
Pctg_array=0.1:0.05:0.96;

fprintf('Iter: \t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\n','Nuclear','|Trace|','N\_Flg','Tol','Res','T\_Flg');

for iter_m=1:size(Pctg_array,2)
    
    N_NumSucc=0;
    T_NumSucc=0;
    
    Pctg=Pctg_array(iter_m);
    for iter_c=1:MaxTest
        
        % signal generation
        f=randsample([0:1/n:(n-1)/n],s-1);
        if f(s-2)+delta<1
            f(s)=f(s-2)+delta;
        else
            f(s)=f(s-2)-delta;
        end
        
        c=rand(s,1)+sqrt(-1)*rand(s,1);
        
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
        
        if trabs<NN
            N_Flg='Succeed';
            N_NumSucc=N_NumSucc+1;
        else
            N_Flg='Fail';
        end
        
        % sampling
        Omega=randsample([1:1:nc-1,nc+1:1:n],round(Pctg*n));

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

        if res_err<tol
            T_Flg='Succeed';
            T_NumSucc=T_NumSucc+1;
        else
            T_Flg='Fail';
        end

        % intermediate data
        fprintf('%4d: \t%.2e\t%.2e\t%8s\t%.2e\t%.2e\t%8s\n',iter_c,NN,trabs,N_Flg,...
            tol,res_err,T_Flg);
    end
    
    N_RateSucc(iter_m)=N_NumSucc/MaxTest;
    T_RateSucc(iter_m)=T_NumSucc/MaxTest;
    
end

% graphical illustration: fix # of sinusoids, illustrate success rate VS # of measurements
figure
plot(Pctg_array,N_RateSucc,'-*',Pctg_array,T_RateSucc,'-o');
xlabel('Percentage of samples');
ylabel('Success rate');
legend('Null space criterion','Tolerance criterion');
