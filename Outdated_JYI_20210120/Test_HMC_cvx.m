clear all
clc

%% save
count=0;    
% save data to xls file
clk=clock;
datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5))]);
filename=strcat('HMC-',datestmp,'-Succ');
filename=strcat(filename,'.xlsx');
col_name={'nc','nSuccess'};
xlswrite(filename,col_name);

%% Problem data
for nc=6:2:6
    nc
    n=2*nc-1;
    k=nc-2;
   
    
    %% Solve EMaC
    trial=1;
    nSucc=0;
    for ii=1:trial
        ii
        %% 1D frequency signal
        %f0=rand(k,1);
        f0=[0.1;0.11;0.12;0.13]
        %f0=[0.1;0.11;0.12]
        x0=zeros(n,1);
        V0=exp(1i*2*pi*kron((0:n-1)',f0'));
        c0=10*(rand(k,1)+1i*rand(k,1));
        x0=V0*c0;
        
        % p=randperm(n);
        Omega=[(1:nc-1)';(nc+1:n)'];
        
        tic
        cvx_solver sdpt3
        cvx_begin sdp quiet
        cvx_precision best
        variable Y(nc,nc) hermitian
        variable Z(nc,nc) hermitian
        variable u(n,1) complex
        Q1=[Y, hankel(u(1:nc),u(nc:end)); hankel(u(1:nc),u(nc:end))', Z];
        minimize 0.5*trace(Y)+0.5*trace(Z)
        subject to
        Q1 >= 0,
        u(Omega) == x0(Omega);
        cvx_end
        toc
        % estimate the frequencies using esprit
        res_err = norm(u-x0)/norm(x0)
        abs_err=norm(u-x0)
        if res_err <= 10^-3
            nSucc=nSucc+1;
        end
    end
    resbuf=[nc,nSucc];
    %% save data to xls file
    count=count+1;
%     xlRange = sprintf('A%d',count);
%     xlsInsert = resbuf;
%     xlswrite(filename,xlsInsert,'Sheet1',xlRange);
end

