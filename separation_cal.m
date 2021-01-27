% This file is to calculate the minimal separation between two atoms, and
% there are only two atoms in the atomic set.
% 
% Based on [Tang, 2015], [Xu, 2017]
%
% Created by JYI, 11/04/2019
% 
%% 
% rng(0)
% N = 127; n = N-1; t = 0:n;
% f = rand(3,1);
% AtomSet = exp(sqrt(-1)*2*pi*kron(t',f'));
% aDist = norm(AtomSet(:,1) - AtomSet(:,2),2);
% SVal = svd(AtomSet);
% aDistD = 2* SVal(end) / sqrt(3);
% 

%%
rng(0); 
N = 2095; n = N-1; t = 0:n; 
figname = sprintf('Atom_separation_%s',datestr(now,'yyyymmddHHMM'));

aNumA = 2;
aDistT = [];
aDistX = [];

for aNum=aNumA
    
    % Tang
    f1 = rand(aNum-1,1); 
    fDist = 2/(2*pi*n); 
    f2 = f1(1)+fDist;
    while f2>1 || f2 <0
        f2 = f1-fDist;
    end

    AtomSet = exp(sqrt(-1)*2*pi*kron(t',[f1;f2]'));

    a1 = AtomSet(:,1);
    a2 = AtomSet(:,end);
    aDist = norm(a1-a2,2); 
    aDistT = [aDistT,aDist];
    fprintf('Tang: Atom separation: %8.4e, Fre 1 %8.4e, Fre 2 %8.4e, Freq separation %8.4e\n', ...
        aDist,f1,f2,fDist);

    % Xu
    fDist = 2/(n);
    f2 = f1(1)+fDist;
    while f2>1 || f2 <0
        f2 = f1-fDist;
    end
    
    AtomSet = exp(sqrt(-1)*2*pi*kron(t',[f1;f2]'));
    Q = rank(AtomSet);
    if Q==aNum
        SVal = svd(AtomSet); % A = [a1,a2]; 
        aDistD = 2* SVal(end) / sqrt(aNum);
        aDistX = [aDistX,aDistD];
    end

    fprintf('Xu: Atom separation %8.4e, Fre 1 %8.4e, Fre 2 %8.4e\n',aDistD,f1,f2);
end

figure; hold on;
plot(aNumA,aDistT,'-*'); plot(aNumA,aDistX,'-o');
legend('[Tang, 2015]','Ours');
xlabel('Number of atoms'); ylabel('Atom separation');

% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,figname,'-dpdf');
% close all;

