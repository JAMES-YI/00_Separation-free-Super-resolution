% This file is to calculate the relation between the frequency distance and
% atom distance for two atoms.
%
% Created by JYI, 03/21/2020
% Modified by JYI, 01/26/2021
% - ToDo: documentation
clear all; clc
%% N = 127
% - use interpolation to estimate the relation between frequency separation
%   and atom separation
% - only consider two atoms

f1 = 0; 
N = 127; t = 0:N-1; dfV = [0.4:-0.001:0.008,0.0079:-0.001:0.003,0.0025:-0.0001:0.0001];
aNumA = 2;
aDistX = [];

% Interpolation between frequency separation and atom separation
for df = dfV
    
    f2 = f1+df;
    AtomSet = exp(sqrt(-1)*2*pi*kron(t',[f1,f2]));
    aDistD = norm(AtomSet(:,1)-AtomSet(:,2));
    aDistX = [aDistX,aDistD];
%    Q = rank(AtomSet);
%     if Q == aNumA
%         SVal = svd(AtomSet); % A = [a1,a2]; 
%         aDistD = 2* SVal(end) / sqrt(aNumA);
%         aDistX = [aDistX,aDistD];
%     end
    
end

%Ind = find(dfV==0.0079); IndT = find(dfV==0.0025); % N=127
figure; hold on;
p1 = plot(dfV,aDistX,'-'); %plot(dfV(Ind),aDistX(Ind),'r*'); plot(dfV(IndT),aDistX(IndT),'b*');
xlabel('$|f_1-f_2|$','Interpreter','latex'); ylabel('$||a(f_1)-a(f_2)||$','Interpreter','latex');

% box on; figname = sprintf('Atom_Dist_VS_Fre_Dist_%d_%s',N,datestr(now,'yyyymmddHHMM'));
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,figname,'-dpdf');
% close all;

% Atom separation bound from frequency separation bound [Tang et al., 2015]
f = [0,0.0004]; 
AtomSet = exp(sqrt(-1)*2*pi*kron(t',[f1,f2]));
aDistD = norm(AtomSet(:,1)-AtomSet(:,2))

%% N = 63

f1 = 0; 
N = 63; t = 0:N-1; dfV = [0.4:-0.001:0.016,0.0159:-0.001:0.006,0.0051:-0.0001:0.0001];
aNumA = 2;
aDistX = [];

% Interpolation between frequency separation and atom separation
for df = dfV
    
    f2 = f1+df;
    AtomSet = exp(sqrt(-1)*2*pi*kron(t',[f1,f2]));
    aDistD = norm(AtomSet(:,1)-AtomSet(:,2));
    aDistX = [aDistX,aDistD];
%     Q = rank(AtomSet);
%     if Q == aNumA
%         SVal = svd(AtomSet); % A = [a1,a2]; 
%         aDistD = 2* SVal(end) / sqrt(aNumA);
%         aDistX = [aDistX,aDistD];
%     end
    
end

%Ind = find(dfV==0.0159); IndT = find(dfV==0.0051); % N=63
p2 = plot(dfV,aDistX,'-'); %plot(dfV(Ind),aDistX(Ind),'r*'); plot(dfV(IndT),aDistX(IndT),'b*');
xlabel('$|f_1-f_2|$','Interpreter','latex'); ylabel('$||a(f_1)-a(f_2)||$','Interpreter','latex');

% Atom separation bound from frequency separation bound [Tang et al., 2015]
f = [0,0.0011]; 
AtomSet = exp(sqrt(-1)*2*pi*kron(t',[f1,f2]));
aDistD = norm(AtomSet(:,1)-AtomSet(:,2))

%% N=3

f1 = 0; 
N = 3; t = 0:N-1; dfV = [0.4:-0.001:0.016,0.0159:-0.001:0.006,0.0051:-0.0001:0.0001];
aNumA = 2;
aDistX = [];

% Interpolation between frequency separation and atom separation
for df = dfV
    
    f2 = f1+df;
    AtomSet = exp(sqrt(-1)*2*pi*kron(t',[f1,f2]));
    aDistD = norm(AtomSet(:,1)-AtomSet(:,2));
    aDistX = [aDistX,aDistD];
%     Q = rank(AtomSet);
%     if Q == aNumA
%         SVal = svd(AtomSet); % A = [a1,a2]; 
%         aDistD = 2* SVal(end) / sqrt(aNumA);
%         aDistX = [aDistX,aDistD];
%     end
    
end

%Ind = find(dfV==0.0159); IndT = find(dfV==0.0051); % N=63
p3 = plot(dfV,aDistX,'-'); %plot(dfV(Ind),aDistX(Ind),'r*'); plot(dfV(IndT),aDistX(IndT),'b*');
xlabel('$|f_1-f_2|$','Interpreter','latex'); ylabel('$||a(f_1)-a(f_2)||$','Interpreter','latex');
legend([p1,p2,p3],'2N-1=127','2N-1=63','2N-1=3');

% Atom separation bound from [Tang et al., 2015]
f = [0,0.167]; 
AtomSet = exp(sqrt(-1)*2*pi*kron(t',[f1,f2]));
aDistD = norm(AtomSet(:,1)-AtomSet(:,2))
%%
% box on; figname = sprintf('Atom_Dist_VS_Fre_Dist_%s',datestr(now,'yyyymmddHHMM'));
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3),fig_pos(4)];
% print(fig,figname,'-dpdf');
% close all;



