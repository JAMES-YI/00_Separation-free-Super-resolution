function [ EDFreq ] = RANM_dual(YonOmega, Omega, N, k, fplot, freq, ck, MaxItr)
%% Reweighted Atoimc Norm Minimization (RANM)
% (NOTE) The frame of this code is obtained from Z. Yang. (https://sites.google.com/site/zaiyang0248/) 
% For displaying the results, few more codes are added without modifying the main code.
%
%
% [Y,u,freq,amp] = ANM_sdpt3(YonOmega, Omega, N, eta)
% 
% ANM_sdpt3 implements the atomic norm minimization problem:
% min_Y ||Y||_A, subject to ||Y_Omega - Yo_Omega||_F <= epsilon
% via duality using SDPT3.
% 
% The dual problem is
% min_{V,H} <V_Omega, Yo_Omega>_R + epsilon*||V_Omega||_F,
% subject to [I V'; V H] >= 0, V_Omegac = 0, and T^*(H) = [1,0,...,0]^T.
% 
% Input:
%   YonOmega: observed measurements on Omega
%   Omega: index set of the measurements
%   N: length of sinusoidal signal of interest
% Output:
%   Y: recovered Y
%   u: u composing the Toeplitz matrix
%   freq: recovered frequency
%   amp: recovered amplitude


Omegac = (1:N)';
Omegac(Omega) = [];

uOld = zeros(N,1);
YOld = zeros(N,1);

for itr = 1:MaxItr
tRANM=tic;
if itr <= 10
    epsilon = 1/2^itr;
else
    epsilon = 1/2^10;
end
W = inv(toeplitz(uOld) + epsilon*eye(N));

% solve the dual problem
cvx_quiet true
cvx_precision default
cvx_solver sdpt3

cvx_begin sdp

  variable X(N+1,N+1) hermitian;
  dual variable U;
  
  X >= 0 : U;
  X(1+Omegac,1) == 0;
  for j = 0:N-1
    sum(diag(W - X(2:N+1,2:N+1), j)) == 0;
  end
  X(1,1) == 1;

  minimize real( trace( YonOmega(Omega)'*X(Omega+1,1) ) )
  
cvx_end
% Y estimate
Y = U(2:1+N,1) * 2;

% postprocessing
u = U(2,2:1+N).';
ED=norm(Y-YOld,2);
if ED < 10^-6
    break;
end
YOld = Y;
uOld = u;

fprintf('RANM runtime/iteration: %f (dis:%f)\n',toc(tRANM),ED);
end

[fEst, amp] = VanDec(u);

amp = amp * 2;

% display Mean Square Error (MSE)
nfEst = max(size(fEst));
fErr = 0;
if nfEst >= k
    freqC = [freq;min(freq)+1;max(freq)-1];
    for i=1:nfEst
        % freq error
        [res, ind] = sort(abs(freqC - fEst(i)),'ascend');
        fErr = fErr + res(1)^2;
    end
else 
    fEstC = [fEst;min(freq)+1;max(freq)-1];
    for i=1:k    
        [res, ~] = sort(abs(freq(i) - fEstC),'ascend');
        fErr = fErr + res(1)^2;
    end
end
EDFreq = sqrt(fErr);
fprintf('RANM runtime/iteration: %f (dis:%f)\n',toc(tRANM),ED);
fprintf('RANM: Euclidean Distance of freq = %f (itr = %d)\n', EDFreq, itr);
fprintf('---------------------------------------------------------\n');

end