% this file is to illustrate performance of identifying underlying frequencies of exponentials after 
% successfully recovering all time samples of signal from its partial observation set
% 
% the signal with has s frequency compoents which form orthogonal atoms, the Hankel matrix has
% fixed dimension N by N
% 
% we will randomly remove one of the atoms on the grid, and instead we add
% another atom which may not be on the grid but is very close to one the
% remainings with distance Delta f=10^(-4)
%
% this file will test 20 times for fixed N to get the successful recovery rate
% 
% the performance of frequencies identification is evaluated by the
% overlapping of underlying frequencies and the local maximizer of imaging
% function
%
% J. Yi on 09/04/2017, jirong-yi@uiowa.edu

% n: # of times instance considered 
% t: time instances considered
% nc: dimension of Hankel matrix
% s: # of exponentials in signal
% delta: smallest frequency separation
% Q: element in null space of sampling operator
% MaxTest: maximal test for each number of measurements
% tol: criterion for determing success of recovery
% mag_min: minimal magnitude of all coefficients

clear all
close all
clc 
warning off

%% Parameters setting up
nc=64;
n=2*nc-1;
t=0:(n-1);
Omega=[1:63,65:n];

a=1;
Q=a*fliplr(diag(ones(nc,1)));
Delta=10^(-4);

mag_min=1;
mag_cl=mag_min/10^2;
mag_rm=1+10.^(0.5*rand(1));
mag_rel=(mag_cl+mag_rm)/mag_min;
s=3;
% 5, 3, 
% s=round((nc^2-2*sqrt(2)*nc*sqrt(nc)-nc*log(nc)+mag_rel*sqrt(2*nc)*log(nc))/(2*(nc+log(nc))));

tol=10^(-3);

fprintf('\t%8s\t%8s\t%8s\t%8s\n','Nuclear','|Trace|','Tol','Res');
        
% mutually orthogonal atoms generation
f=randsample([0:1/nc:(nc-1)/nc],s);
c_mag=[mag_min;mag_rm;1+10.^(0.5*rand(s-2,1))];
c_pha=rand(s,1)*2*pi;

% remove and add
Ind=randsample(s,2);
Rm=Ind(1);
Cl=Ind(2);

if (f(Cl)+Delta>=1)
    f_cl=f(Cl)-Delta;
else
    f_cl=f(Cl)+Delta;
end

c_cl_mag=mag_cl;
c_cl_pha=rand(1)*2*pi;

f_p=[f(1:Rm-1),f(Rm+1:end),f_cl];
c_p_mag=[c_mag(1:Rm-1);c_mag(Rm+1:end);c_cl_mag];
c_p_pha=[c_pha(1:Rm-1);c_pha(Rm+1:end);c_cl_pha];
c_p=c_p_mag.*exp(sqrt(-1)*c_p_pha);

% perturbed signal 
V=exp(1i*2*pi*kron(t',f_p));
x=V*c_p; 

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

% intermediate data
fprintf('\t%.2e\t%.2e\t%.2e\t%.2e\n',NN,trabs,tol,res_err);

if res_err<tol && trabs<NN
    
    fprintf('All time domain samples have been recovered successfully!\n');
    fprintf('Underlying frequencies are being identified...\n');
    
    % images for all frequencies in [0,1)
    hank_u=hankel(u(1:nc),u(nc:n));
    [U_u,D_u,V_u]=svd(hank_u);
    U2_u=U_u(:,s+1:end);
    w=0:0.00001:1;
    
    J_w=zeros(1,length(w));
    for iter_w=1:length(w)
        
        Atm_w=exp(sqrt(-1)*2*pi*w(iter_w)*[0:nc-1]');
        J_w(iter_w)=(norm(Atm_w,2))/(norm(U2_u'*Atm_w,2));
        
    end
    
    % images for underlying frequencies of exponentials
    J_f=zeros(1,length(f_p));
    for iter_f=1:length(f_p)
        
        Atm_f=exp(sqrt(-1)*2*pi*f_p(iter_f)*[0:nc-1]');
        J_f(iter_f)=(norm(Atm_f,2))/(norm(U2_u'*Atm_f,2));
        
    end
    
    figure 
    hold on
    plot(w,J_w);
    xlabel('Frequency w');
    ylabel('Imaging J(w)');
    title('Frequencies Identification');
    plot(f_p,zeros(size(f_p)),'o');
    
    
%     subplot(1,2,1)
%     plot(w,J_w);
%     xlabel('Frequency w');
%     ylabel('Imaging J(w)');
%     title('Imaging Over [0,1)');
    
%     subplot(1,2,2)
%     plot(f_p,J_f,'o');
%     xlabel('Frequency w');
%     ylabel('Imaging J(w)');
%     title('Imaging Over Underlying Frequencies');
    
    fprintf('Frequency identification is done!\n');
    
end