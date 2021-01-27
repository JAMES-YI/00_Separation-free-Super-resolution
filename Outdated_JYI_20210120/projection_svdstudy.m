n=6;
h=0.5;
basis=dftmtx(n)
s=3;
v=randn(n, s)+sqrt(-1)*randn(n,s);
A=fliplr(eye(n))
B=v*v';
[u1,C,v1]=svd(B);
u1=(1/sqrt(n))*basis(:,[3,4,6]);
v1=(1/sqrt(n))*basis(:,[3,4,6]);


%b=u1(:,s+1:n)'*A*conj(v1(:,(s+1):n));
b=u1'*A*conj(v1);

%clear all
%clc

% % %% Problem data 
% % nc=3;
% % n=2*nc-1;
% % k=nc-1;
% % m=2*nc-1;
% % 
% % %% 1D frequency signal
% % f0=rand(k,1);
% % x0=zeros(m,1);
% % V0=exp(1i*2*pi*kron((0:n-1)',f0'));
% % %V0=real(exp(1i*2*pi*kron((0:n-1)',f0')));
% % 
% % c0=rand(k,1)+1i*rand(k,1);
% % %c0=real(rand(k,1)+1i*rand(k,1));
% % x0=V0*c0;
% % 
% % p=randperm(n);
% % %Omega=p(1:m)';
% % Omega=[1:nc-1,nc+1:n].';
% % %%%55555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











% A=randn(n, n)
% svd(A);
% B=A(2:n,2:n);
% svd(B);
% 
% r=2;
% A=randn(n,n)+sqrt(-1)*randn(n,n);
% [U,S,V]=svd(A)
% for i=r+1:n
%     S(i,i)=0;
% end
% A=U*S*V';
% Theta=[real(U(:,1:r)), -imag(U(:,1:r));imag(U(:,1:r)), real(U(:,1:r))];
% Gamma=[real(V(:,1:r)), -imag(V(:,1:r));imag(V(:,1:r)), real(V(:,1:r))];
% TNull=null(Theta.');
% GNull=null(Gamma.');
% B=randn(2*(n-r),2*(n-r));
% sv=svd(B);
% norms=sv(1,1);
% B=B/norms;
% Delta=TNull*B*GNull.';
% RDelta=zeros(size(Delta));
% RDelta(1:n,1:n)=(Delta(1:n,1:n)+Delta((n+1):2*n,(n+1):(2*n)))/2;
% RDelta((n+1):2*n, (n+1):2*n)=(Delta(1:n,1:n)+Delta((n+1):2*n,(n+1):(2*n)))/2;
% RDelta((n+1):2*n, 1:n)=(Delta((n+1):2*n,1:n)-Delta(1:n,(n+1):(2*n)))/2;
% RDelta(1:n,(n+1):2*n)=-(Delta((n+1):2*n,1:n)-Delta(1:n,(n+1):(2*n)))/2;
% FirstPart=Theta(1:n,1:r)*Gamma(1:n,1:r).'+Theta((n+1):2*n,1:r)*Gamma((n+1):2*n,1:r).'
% SecondPart=Theta((n+1):2*n,1:r)*Gamma(1:n,1:r).'-Theta(1:n,1:r)*Gamma((n+1):2*n,1:r).'
% Com=FirstPart+sqrt(-1)*SecondPart;
% W=RDelta(1:n, 1:n)+sqrt(-1)*RDelta((n+1):2*n, 1:n);
% U(:,1:r)'*W
% W*V(:,1:r)
% svd(W)
% 

%%%%%%%%%%test program for singular vectors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_col=0;

for trial=1:10 % originally 100000

%% Problem data 
nc=3;
n=2*nc-1;
k=nc-1;
m=2*nc-2;

%% 1D frequency signal
f0=rand(k,1);
x0=zeros(m,1);
V0=exp(1i*2*pi*kron((0:n-1)',f0'));
%V0=real(exp(1i*2*pi*kron((0:n-1)',f0')));

c0=rand(k,1)+1i*rand(k,1);
%c0=real(rand(k,1)+1i*rand(k,1));
x0=V0*c0;

p=randperm(n);
%Omega=p(1:m)';
Omega=[1:nc-1,nc+1:n].';
%%%55555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=hankel(x0(1:nc),x0(nc:end))
[u1,C,v1]=svd(H);
A=fliplr(eye(nc))
b=trace((u1*(v1)')'*A)
s_col=s_col+b;

x1=V0*(c0+0.009* randn(k,1))
H2=hankel(x1(1:nc),x1(nc:end))
[u2,C2,v2]=svd(H2);
A=fliplr(eye(nc))
b2=trace((u2*(v2)')'*A)
%s_col2=s_col2+b2;

end
s_col/100000