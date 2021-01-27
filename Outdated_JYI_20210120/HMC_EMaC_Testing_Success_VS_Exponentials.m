% this file is to illustrate successful recovery rate for signal with two
% frequency compoents close to each other, with frequency separation delta=0.001

% n: number of time instance we consider
% t: time instances we consideration
% nc: dimension of Hankel matrix
% s: # of sinusoid components, the maximal rank the Hankel matrix can
% achieve
% s_array: # of exponentials in the signal, each element in s_array
% represent a case
% delta: frequency separation
% f: frequencies vector
% V: Fourier matrix
% Omega: observed sample indices
% c: coefficients of sinusoid components
% x: true signal 
% u: recovered signal
% tol: criterion for successful recovery
% MaxItr_c: maximal # of choices for coefficients
% MaxItr_r: maximal # of tests for each group of coefficients
% NumSucc: # of successful recovery
% Succ_Rate: successful recovery rate
% a: scaling factor for element from null space of projector
% 

% J. Yi on 05/30/2017, jirong-yi@uiowa.edu


clear all
close all
clc 
warning off

%% Parameters setting up
n=127;
% n=11;
t=0:(n-1);
nc=(n+1)/2;
s_array=3:4:round(0.3*n);
a=1;

delta=1/256;
% delta: 0.1,0.05,0.01,0.001

MaxTest=10;
tol=10^(-6);
Pctg=0.4;

fprintf('Iter: \t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\n','Nuclear','|Trace|','N\_Flg','Tol','Res','T\_Flg');

for iter_m=1:size(s_array,2)
    
    N_NumSucc=0;
    T_NumSucc=0;
   
    Q=a*fliplr(diag(1,-(nc-1)));
    s=s_array(iter_m);
    
    for iter_c=1:MaxTest
        
        f=randsample([0:1/n:(n-1)/n],s-1);
        if f(s-2)+delta<1
            f(s)=f(s-2)+delta;
        else
            f(s)=f(s-2)-delta;
        end
        
        c=rand(s,1)+sqrt(-1)*rand(s,1);

        V=exp(1i*2*pi*kron(t',f));
        % c=100*rand(s,1)+100i*rand(s,1);
        %     cReal=10*rand(s,1);
        %     c(1,1)=cReal(1)+sqrt(-1)*cReal(1)/10;
        %     c(2,1)=cReal(2)+sqrt(-1)*cReal(2)*10;

        x=V*c; 

        Omega=randsample([1:1:n-1],round(Pctg*n));

        % signal recovery by EMaC
%         tic
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
%         toc;

        % results report

        % if res_err<tol
        %     Flg='Succeed';
        %     NumSucc=NumSucc+1;
        % else
        %     Flg='Fail';
        % end
        
        % null space criterion

        hank_x=hankel(x(1:(n+1)/2),x((n+1)/2:n));
        [U_x,D_x,V_x]=svd(hank_x);

    %     U_b=U_x(:,(n+1)/2);
    %     V_b=V_x(:,(n+1)/2);
    %     U_x=U_x(:,1:(n-1)/2);
    %     V_x=V_x(:,1:(n-1)/2);

        U_b=U_x(:,s+1:end);
        V_b=V_x(:,s+1:end);
        U_x=U_x(:,1:s);
        V_x=V_x(:,1:s);

        tr(iter_c)=trace(V_x*U_x'*Q);
        trabs(iter_c)=abs(tr(iter_c));
        NN(iter_c)=sum(svd(U_b'*Q*V_b));
        if trabs(iter_c)<NN(iter_c)
            N_Flg='Succeed';
            N_NumSucc=N_NumSucc+1;
        else
            N_Flg='Fail';
        end

        % tolerance criterion
        res_err(iter_c) = norm(u-x)/norm(x);

        if res_err(iter_c)<tol
            T_Flg='Succeed';
            T_NumSucc=T_NumSucc+1;
        else
            T_Flg='Fail';
        end

        % intermediate data
        fprintf('%4d: \t%.2e\t%.2e\t%8s\t%.2e\t%.2e\t%8s\n',iter_c,NN(iter_c),trabs(iter_c),N_Flg,...
            tol,res_err(iter_c),T_Flg);
    end
    
    N_RateSucc(iter_m)=N_NumSucc/MaxTest;
    T_RateSucc(iter_m)=T_NumSucc/MaxTest;
    
end

% graphical illustration: fix # of sinusoids, illustrate success rate VS # of measurements
figure
plot(s_array,N_RateSucc,'-*',s_array,T_RateSucc,'-o');
xlabel('Number of exponentials');
ylabel('Success rate');
legend('Null space criterion','Tolerance criterion');


    
% RCount=RCount+1;
% xlswrite(filename,{'SuccessRate',RateSucc},sheetname,sprintf('F%d:G%d',RCount,RCount));

% graphical illustration

% figure
% subplot(3,1,1)
% plot(trabs);
% title('|trace|');
% 
% subplot(3,1,2)
% plot(NN);
% title('NN');
% 
% subplot(3,1,3)
% plot(NN-trabs);
% title('NN-|trace|');
% 
% figure
% subplot(3,1,1)
% hist(trabs);
% title('|trace|');
% 
% subplot(3,1,2)
% hist(NN);
% title('NN');
% 
% subplot(3,1,3)
% hist(NN-trabs);
% title('NN-|trace|');

    


%     cRange=sprintf('A%d',RCount);
%     data_c={'c1',real(c(1)),imag(c(1));'c2',real(c(2)),imag(c(2))};
%     xlswrite(filename,data_c,'sheet1',cRange);
%     
%     sRange=sprintf('E%d',RCount);
%     xlswrite(filename,[NumSucc,RateSucc],'sheet1',sRange);
%     
%     RCount=RCount+3;
%     
%     clear c Y Z Q1 u NumSucc 
    




%% Solve EMaC
% this part folves the following optimization
% min_Y,Z,u 0.5(tr(Y)+tr(Z))
% s.t. [Y,hankel(u(1:nc),u(nc:end));hankel(u(1:nc),u(nc:end))',Z] >=0
%      u(Omega)=x_0(Omega)
%      Y(nc,nc), Z(nc,nc) are Hermitian, u(n,1) is complex

% tic
% cvx_solver sdpt3
% cvx_begin sdp 
%     variable Y(nc,nc) hermitian
%     variable Z(nc,nc) hermitian
%     variable u(n,1) complex 
%     Q1=[Y, hankel(u(1:nc),u(nc:end)); hankel(u(1:nc),u(nc:end))', Z];
%     minimize 0.5*trace(Y)+0.5*trace(Z)
%     subject to
%         Q1 >= 0,
%         u(Omega) == x0(Omega);
%     
% cvx_end
% toc
% % estimate the frequencies using esprit
% res_err = norm(u-x0)/norm(x0)

