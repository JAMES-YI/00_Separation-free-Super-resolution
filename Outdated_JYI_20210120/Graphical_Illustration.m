% this file is to illustrate the results of 1000 cases where the
% coefficients are randomly generated and the frequencies are fixed

close all
clc

%% loading trace data
% filename='data_U, V, and Tr(VU^HQ)_26-May-2017-16h17m.xlsx';
% sheetname='U&V&Tr';

% tr=[]
% for i=3:4:3999
%     fprintf('%d\n',i)
%     Range=sprintf('N%d',i);
%     tr=[tr,xlsread(filename,sheetname,Range)];
% end
% 
% data=load('tr.mat');
% tr=data.tr;
% 
% figure
% plot(tr);
% title('tr(VU^HQ)');
% xlabel('Case number');
% ylabel('Trace');

%% probability plot

t=-20:0.5:20;

% fixed a=1+1i and c=1, different number of time sample instances N=5, 50, 500,
% 5000, 50000
a=1+1i;
c=1;

for N=[5,50,500,5000,50000]
    K=N-c*sqrt(log(N)*N);
    Prob=3*exp(-t.^2*0.5./(K+t/3));
    plot(t,Prob);
    hold on;
end

title('Probability of trace');
ylabel('Prob(|tr(VU^HQ)|)\geq |a|t');
xlabel('t');
legend('N=5','N=50','N=500','N=5000','N=50000');


% fixed a=1+1i and N=50000, different c=1,10,20,30,40
figure
clear a c N

a=1+1i;
N=50000;

for c=[1,10,20,30,40]
    K=N-c*sqrt(log(N)*N);
    Prob=3*exp(-t.^2*0.5./(K+t/3));
    plot(t,Prob);
    hold on;
end

title('Probability of trace');
ylabel('Prob(|tr(VU^HQ)|)\geq |a|t');
xlabel('t');
legend('c=1','c=10','c=20','c=30','c=40');

% fixed c=1 and N=50000. different a=1+1i, 10+10i, 20+20i, 30+30i, 40+40i
% no need





