function [ A ] = generationFourierMatrixAtomByAtom( Nsample, Ndelta_f )
%% n x m Fourier matrix generation
n = Nsample; % row
m = Ndelta_f; % column
delta_f = 1/Ndelta_f;
A = zeros(n,m);
atom = @(ii) exp(+1*1i*2*pi*delta_f*(ii-1)*[0:n-1]');
for jj=1:m
    A(:,jj)=atom(jj);
end

end

