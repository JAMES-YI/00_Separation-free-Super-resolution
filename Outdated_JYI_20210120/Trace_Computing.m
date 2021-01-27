% this file is to illustrate the row and column subspaces of Hankel matrix
% generated from a spectrally sparse signal

% for anti-diagonal: from top left to right bottom, labeled from 1st up to
% 2n+1-th anti-diagonal

% s: # of sinusoid
% omega: frequencies
% 2n+1: # of time instances
% t: time instances of signal
% c: coefficients of sinusoids
% U, V: n+1 by s
% tr: absolute value of trace of VU^HQ
% delta: frequency separation
% a: scaling factor in null space of operator P_\Omega
% Omega: sampling set {0,1,2,4,5}

% by J. Yi 05/18/2017, jirong-yi@uiowa.edu

% 

warning off
clc
close all

s=2;
delta=[0.1,0.01,0.001,0.0001];
a=(1+sqrt(-1));
% omega(1)=randn(1);
% omega(2)=omega(1)+1/6;
% omega(1)=randn(1);
% omega(2)=omega(1)+1;
n=2;
t=[0:1:2*n]';
caseN=1000;

% RCount=3;
% clk=clock;
% datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5)),'m']);
% filename = strcat('data_U, V, and Tr(VU^HQ)_',datestmp,'.xlsx');
% matfile=strcat('tr_',datestmp,'.mat');
% 
% sheetname='U&V&Tr';
% columnname={'coefficients','real part','imaginary part',' ','column','space','U',...
%     ' ','row','space','V',' ','missing anti-diagonal','trace |tr(VU^HQ)|'};
% xlswrite(filename,columnname,sheetname,'A1:M1');

omega(1)=rand(1);
for iter0=1:size(delta,2)
    if omega(1)+delta(iter0)<1
        omega(2)=omega(1)+delta(iter0);
    else
        omega(2)=omega(1)-delta(iter0);
    end
    
    for iter=1:caseN
        fprintf('Frequency separation %d, Case %d...',delta(iter0),iter);
        sig_x=zeros(size(t,1),1);
        c=rand(s,1).*exp(sqrt(-1)*2*pi*randn(s,1));
%         c(2)=100*c(2);
        for j=1:1:s
            sig_x=sig_x+c(j)*exp(sqrt(-1)*2*pi*omega(j).*t);
        end

        hank_x=hankel(sig_x(1:n+1),sig_x(n+1:2*n+1));
        [U_x,D_x,V_x]=svd(hank_x);
        
        U_b=U_x(:,n+1);
        V_b=V_x(:,n+1);
        U_x=U_x(:,1:n);
        V_x=V_x(:,1:n);

        Q=a*fliplr(eye(3));
        tr(iter)=trace(V_x*U_x'*Q);
        trabs(iter)=abs(trace(V_x*U_x'*Q));
        NN(iter)=sum(svd(U_b'*Q*V_b));
        
        Diff(iter)=NN(iter)-trabs(iter);

    %     columns=[1:1:2*n+1];
    %     a=1+1i;
    %     tr=[];
    %     for d=columns
    %         if d<=n+1
    %             Q1=a*diag(ones(d,1),n+1-d);
    %         else
    %             Q1=a*diag(ones(2*n+2-d,1),n+1-d);
    %         end
    %         Q=fliplr(Q1);
    %         tr=[tr,abs(trace(V_x*U_x'*Q))];
    %     end

    %     cRange=sprintf('A%d',RCount);
    %     data_c={'c1',real(c(1)),imag(c(1));'c2',real(c(2)),imag(c(2))};
    %     xlswrite(filename,data_c,sheetname,cRange);
    %     
    %     uRange=sprintf('E%d',RCount);
    %     xlswrite(filename,U_x,sheetname,uRange);
    %     
    %     vRange=sprintf('I%d',RCount);
    %     xlswrite(filename,V_x,sheetname,vRange);
    %     
    %     trRange=sprintf('M%d',RCount);
    %     xlswrite(filename,{'3-rd',tr(iter)},sheetname,trRange);
    %     
    %     RCount=RCount+4;

        fprintf('done!\n');
    end
    
    figure 
    subplot(6,1,1)
    plot([1:1:caseN]',trabs,'-');
    title('|trace|');
    
    subplot(6,1,2)
    plot([1:1:caseN]',NN,'-');
    title('NN');
    
    subplot(6,1,3)
    plot([1:1:caseN]',Diff,'-');
    title('NN-|trace|');
    
    subplot(6,1,4)
    hist(trabs);
    title('|trace|');
    
    subplot(6,1,5)
    hist(NN);
    title('NN');
    
    subplot(6,1,6)
    hist(Diff);
    title('NN-|trace|');
    

end
% 
% figure
% plot(tr);
% xlabel('Coefficients case');
% ylabel('Tr(VU^HQ)');
% title('|c_1|<<|c_2|');





