% this file is to illustrate successful recovery rate for signal with two
% frequency compoents close to each other, with frequency separation delta=0.001

% n: number of time instance we consider
% t: time instances we consideration
% nc: dimension of Hankel matrix
% s: # of sinusoid components, the maximal rank the Hankel matrix can
% achieve
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

%% Parameters setting up
n=5;
% n=11;
t=0:(n-1);
nc=(n+1)/2;
s=2;
a=1;

delta=0.001;
% delta: 0.1,0.05,0.01,0.001
tol=10^(-10);
% tol=10^-4, 10^-6, 10^-8, 10^-10

Omega=[1,2,4,5]'; % 100% successful rate
% Omega=[1,2,5]'; 100% successful rate
% Omega=[1,5]'; 100% success rate
% Omega=[1,3,5,8,10]';

% clk=clock;
% datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5)),'m']);
% filename = strcat('HMC by EMaC',datestmp,'.xlsx');
% columnname={'Coefficients','Real part','Imaginary part',' ','Success number','Success rate'};
% xlswrite(filename,columnname,'sheet1','A1:F1');

% filename=strcat('Separation and Precision');
% sheetname=strcat('Separa_',num2str(delta),',',' Precs_',num2str(tol));
% columnname={'Test','NN','|trace|','trReal','trImageinary','RelaResidual','Status'};
% xlswrite(filename,columnname,sheetname,'A1:G1');

filename='Recovery Performance (Null space condition)';
sheetname=strcat('Separa_',num2str(delta),',',' Precs_',num2str(tol));
columnname={'Test','NN','|trace|','trReal','trImageinary','RelaResidual','Status'};
xlswrite(filename,columnname,sheetname,'A1:G1');

RCount=3;
MaxItr_c=100;
NumSucc=0;

for iter_c=1:MaxItr_c

f(1)=rand(1);
if f(1)+delta>1
    f(2)=f(1)-delta;
else
    f(2)=f(1)+delta;
end

V=exp(1i*2*pi*kron(t',f));
c=100*rand(s,1)+100i*rand(s,1);
x=V*c;

% signal recovery by EMaC

fprintf('Test # %d\n',iter_c);

tic
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
toc;

% results report

res_err(iter_c) = norm(u-x)/norm(x);

% if res_err<tol
%     Flg='Succeed';
%     NumSucc=NumSucc+1;
% else
%     Flg='Fail';
% end

hank_x=hankel(x(1:(n+1)/2),x((n+1)/2:n));
[U_x,D_x,V_x]=svd(hank_x);

U_b=U_x(:,(n+1)/2);
V_b=V_x(:,(n+1)/2);
U_x=U_x(:,1:(n-1)/2);
V_x=V_x(:,1:(n-1)/2);

Q=a*fliplr(eye(3));
tr(iter_c)=trace(V_x*U_x'*Q);
trabs(iter_c)=abs(tr(iter_c));
NN(iter_c)=sum(svd(U_b'*Q*V_b));

if trabs(iter_c)<NN(iter_c)
    Flg='Succeed';
    NumSucc=NumSucc+1;
else
    Flg='Fail';
end

% saving data

% columnname={'Test','NN','|trace|','trReal','trImageinary','RelaResidual','Status'};
columndata={iter_c,NN(iter_c),trabs(iter_c),imag(tr(iter_c)),real(tr(iter_c)),res_err(iter_c),Flg};
xlswrite(filename,columndata,sheetname,sprintf('A%d:G%d',RCount,RCount));

RCount=RCount+1;

end
    
RateSucc=NumSucc/MaxItr_c;

RCount=RCount+1;
xlswrite(filename,{'SuccessRate',RateSucc},sheetname,sprintf('F%d:G%d',RCount,RCount));

% graphical illustration

figure
subplot(3,1,1)
plot(trabs);
title('|trace|');

subplot(3,1,2)
plot(NN);
title('NN');

subplot(3,1,3)
plot(NN-trabs);
title('NN-|trace|');

figure
subplot(3,1,1)
hist(trabs);
title('|trace|');

subplot(3,1,2)
hist(NN);
title('NN');

subplot(3,1,3)
hist(NN-trabs);
title('NN-|trace|');

    

%     
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

