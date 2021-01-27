for i=1:1000000

Delta=randn(2*n,2*n);
a=svd(Delta);
RDelta=Delta;%zeros(size(Delta));
%RDelta(1:n,1:n)=(Delta(1:n,1:n)+Delta((n+1):2*n,(n+1):(2*n)))/2;
%RDelta((n+1):2*n, (n+1):2*n)=(Delta(1:n,1:n)+Delta((n+1):2*n,(n+1):(2*n)))/2;
RDelta((n+1):2*n, 1:n)=(Delta((n+1):2*n,1:n)-Delta(1:n,(n+1):(2*n)))/2;
RDelta(1:n,(n+1):2*n)=-(Delta((n+1):2*n,1:n)-Delta(1:n,(n+1):(2*n)))/2;
b=svd(RDelta);
if(b(1,1)>a(1,1))
   111111111111111111111111111111111111111111
end

end