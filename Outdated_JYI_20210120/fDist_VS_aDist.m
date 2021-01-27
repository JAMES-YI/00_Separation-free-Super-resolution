% This file is to empirically calculate the relation between atom distance
% and frequency distance.
% 
% Created by JYI, 03/30/2020
% - 
% 
%%
% - time samples: consider only samples at 0,32,64,93,127
% - atom frequencies for calculating atom distance using [Xu et al., 2017]: 0/N, 32/N, 64/N, 93/N, 127/N
% - atom frequencies for calculating atom distance using [Tang et al.,
%   2015]: 0, 2/((N-1)*2*pi)
% - 

rng(0); 
aDistT = []; aDistX = [];
NVec = 127; 

for N=NVec
    n = N-1; t = [0,32,64,93,127]; 
    figname = sprintf('Atom_separation_%s',datestr(now,'yyyymmddHHMM'));

    f = [1,3,5,7,9]./N; aNumA = numel(f); 
    AtomSet = exp(sqrt(-1)*2*pi*kron(t',f));

    % Tang
    fTang = [0,2/(n*2*pi)]+0.4;
    AtomSetTang = exp(sqrt(-1)*2*pi*kron(t',fTang));
    aDist = norm(AtomSetTang(:,1)-AtomSetTang(:,2),2);
    aDistT = [aDistT,aDist];
    fprintf('Tang: Atom separation %4.2f\n',aDist);
    % Xu
    Q = rank(AtomSet);
    if Q==aNumA
        SVal = svd(AtomSet); % A = [a1,a2]; 
        aDistD = 2* SVal(end) / sqrt(aNumA);
        aDistX = [aDistX,aDistD];
    else
        SVal = svd(AtomSet(:,1:Q)); % take the first Q atoms
        aDistD = 2*SVal(end) / sqrt(Q);
        aDistX = [aDistX,aDistD];
    end

    fprintf('Xu: Atom separation %8.4e\n',aDistD);
end

%% 
% - two atoms in total
% - atom frequencies: 0, df
% - atoms consist of samples are specific locations: 0,32,64,93,127
% - give plot: atom distance VS frequency distance

f1 = 0; 
N = 127; t = [0,32,64,93,127]; dfV = [0.4:-0.0001:0.0001];
aNumA = 2; aDistX = [];

for df = dfV
    
    f2 = f1 + df;
    AtomSet = exp(sqrt(-1)*2*pi*kron(t',[f1,f2]));
    aDistD = norm(AtomSet(:,1)-AtomSet(:,2));
    aDistX = [aDistX,aDistD];
    
end

%Ind = find(dfV==0.0079); IndT = find(dfV==0.0025); % N=127
figure; hold on;
p1 = plot(dfV,aDistX,'-'); %plot(dfV(Ind),aDistX(Ind),'r*'); plot(dfV(IndT),aDistX(IndT),'b*');
xlabel('$|f_1-f_2|$','Interpreter','latex'); ylabel('$||a(f_1)-a(f_2)||$','Interpreter','latex');