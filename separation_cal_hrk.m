% This file is to calculate the minimal separation between two atoms, and
% there are only two atoms in the atomic set.
% 
% Based on [Tang, 2015], [Xu, 2017]
%
% Created by JYI, 11/15/2019
% 

%%
rng(0); 
aDistT = []; aDistX = [];
NVec = 3; % 255, 127, 63, 31, 15, 7; [255,127,63,31,15,7,3]

for N=NVec
    n = N-1; t = 0:n; 
    figname = sprintf('Atom_separation_%s',datestr(now,'yyyymmddHHMM'));

    f = t/N; fDist = 1/N; aNumA = numel(f); 
    AtomSet = exp(sqrt(-1)*2*pi*kron(t',f));

    % Tang

%     a1 = AtomSet(:,1);
%     a2 = AtomSet(:,2);
%     aDist = norm(a1-a2,2); 
%     aDistT = [aDistT,aDist];
%     fprintf('Tang: Atom separation: %8.4e, Freq separation %8.4e\n',aDist,fDist);

%     n = NVec-1; t=0:n;
    fTang = [0,2/(n*2*pi)]+0.9;
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
    end

    fprintf('Xu: Atom separation %8.4e\n',aDistD);
end

figure; hold on;
plot(NVec,aDistT,'-*'); 
plot(NVec,aDistX,'-o');
legend('[Tang, 2015]','Ours');
xlabel('2N-1'); ylabel('Atom separation');
box on;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3),fig_pos(4)];
print(fig,figname,'-dpdf');
close all;

