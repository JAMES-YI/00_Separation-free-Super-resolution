n=6;
% h=0.5;
% v=h.^([1:1:n])
% A=fliplr(eye(n))
% B=v'*v;
% [u1,C,v1]=svd(B);
% u1(:, 2:n)'*A*v1(:, 2:n);
% 
% 
% A=randn(n, n)
% svd(A);
% B=A(2:n,2:n);
% svd(B);

r=2;
A=randn(n,n)+sqrt(-1)*randn(n,n);
[U,S,V]=svd(A)
for i=r+1:n
    S(i,i)=0;
end
A=U*S*V';
Theta=[real(U(:,1:r)), -imag(U(:,1:r));imag(U(:,1:r)), real(U(:,1:r))];
Gamma=[real(V(:,1:r)), -imag(V(:,1:r));imag(V(:,1:r)), real(V(:,1:r))];
TNull=null(Theta.');
GNull=null(Gamma.');
B=randn(2*(n-r),2*(n-r));
sv=svd(B);
norms=sv(1,1);
B=B/norms;
Delta=TNull*B*GNull.';
RDelta=zeros(size(Delta));
RDelta(1:n,1:n)=(Delta(1:n,1:n)+Delta((n+1):2*n,(n+1):(2*n)))/2;
RDelta((n+1):2*n, (n+1):2*n)=(Delta(1:n,1:n)+Delta((n+1):2*n,(n+1):(2*n)))/2;
RDelta((n+1):2*n, 1:n)=(Delta((n+1):2*n,1:n)-Delta(1:n,(n+1):(2*n)))/2;
RDelta(1:n,(n+1):2*n)=-(Delta((n+1):2*n,1:n)-Delta(1:n,(n+1):(2*n)))/2;
FirstPart=Theta(1:n,1:r)*Gamma(1:n,1:r).'+Theta((n+1):2*n,1:r)*Gamma((n+1):2*n,1:r).'
SecondPart=Theta((n+1):2*n,1:r)*Gamma(1:n,1:r).'-Theta(1:n,1:r)*Gamma((n+1):2*n,1:r).'
Com=FirstPart+sqrt(-1)*SecondPart;
W=RDelta(1:n, 1:n)+sqrt(-1)*RDelta((n+1):2*n, 1:n);
U(:,1:r)'*W
W*V(:,1:r)
svd(W)

