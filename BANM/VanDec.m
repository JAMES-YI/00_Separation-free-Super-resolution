function [freq, amp, sigma_in_u] = VanDec(u)
%
% [freq, amp, sigma_in_u] = VanDec(u)
% VanDec implements the Vandermonde decomposition given the solution u.
% The result is sorted in the descending order of amp
% 
% References: 
%   Blu et al, "Sparse sampling of signal innovations," IEEE Signal Processing Magazine, 2008
%   Z. Yang and L. Xie, "On gridless sparse methods for line spectral estimation from ...
%      complete and incomplete data," http://arxiv.org/abs/1407.2490, May 2014.
% 
% Written by Z. Yang, Aug 2013

prec = 1e-4;

toepu = (toeplitz(u)+toeplitz(u)')/ 2;
eigval = sort(real(eig(toepu)));
if eigval(1) > 0
    sigma_in_u = eigval(1);
    u(1) = u(1) - sigma_in_u;
    eigval = eigval - sigma_in_u;
else
    sigma_in_u = 0;
end

K = sum(eigval > prec*eigval(end));

M = length(u);

% solve h
h = toeplitz([conj(u(2)); u(1:M-1)], conj(u(2:K+1))) \ (-u(1:M));

% solve the zeros of H(z) in the form of exp(-1i * 2 * pi * theta)
r = roots([1; h]);
r = r ./ abs(r);

% solve the amper vector
mat = flipud(vander([r; zeros(M-K,1)]).');
% amp = mat(:, 1:K) \ u;
amp = [real(mat(:, 1:K)); imag(mat(:, 1:K))] \ [real(u); imag(u)];
[amp, idx] = sort(real(amp), 'descend');
Khat = sum(amp>0);
amp = amp(1:Khat);

% determine the frequency
freq = - phase(r(idx(1:Khat))) / (2*pi);
freq = mod(freq, 1);

end