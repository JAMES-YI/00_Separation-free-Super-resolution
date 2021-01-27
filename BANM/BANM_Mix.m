function [ EDFreq ] = BANM_Mix( xstar, J, N, k, n, fplot, freq, ck, blkRange, delta_f )
%% Block Iterative Reweighted L1 and Atomic Norm Minimization Mixture (BANM-Mix) 
% (NOTE) Atomic Norm Minimization with block prior code is used. 
%
% For detailed information, please see our previous paper:
% Spectral Super-resolution with prior information by K. V. Mishra, et al.
%
% Myung (Michael) Cho
% myung-cho@uiowa.edu

%% Changing index-wise block to continuous block in [0,1] 
tBANM=tic;
alpha_norm = [];
beta_norm = [];
for ii=1:N
    if ii == 1
        if blkRange(ii) == 1 && blkRange(N) == 0
            alpha_norm = [alpha_norm;1/N];
        end
    else
        if (blkRange(ii-1) == 0) && (blkRange(ii) == 1)
            if ii == N && blkRange(1) == 0
                beta_norm = [beta_norm;ii/N];
            else
                alpha_norm = [alpha_norm;ii/N];
            end
        end
        if ii == N
            if (blkRange(ii-1) == 1) && (blkRange(ii) == 1) && (blkRange(1) == 0)
                beta_norm = [beta_norm;ii/N];
            end
        end
        if (blkRange(ii-1) == 1) && (blkRange(ii) == 0)
            beta_norm = [beta_norm;ii/N];
        end
    end
end
if ~(isempty(beta_norm)) && ~(isempty(alpha_norm)) 
if (beta_norm(1) < alpha_norm(1))
    beta_norm = [beta_norm(2:end);beta_norm(1)+1];
end


alpha = alpha_norm.*2*pi;
beta  = beta_norm.*2*pi;

% Consider tau situation (tau value)
% For detailed information about tau value, 
% please see the equation (III.28) of the following paper: 
% Spectral Super-resolution with prior information by K. V. Mishra, et al.
p = length(alpha); % Number of frequency bands
for ii=1:p
    if (alpha_norm(ii) < 0.5) && (beta_norm(ii) > 0.5)
        tau(ii) = 2*pi*(beta_norm(ii) - 0.5 + 0.001);
        alpha(ii) = alpha(ii) - tau(ii);
        beta(ii) = beta(ii) - tau(ii);
    elseif (alpha_norm(ii) < 1.5) && (beta_norm(ii) > 1.5)
        tau(ii) = 2*pi*(beta_norm(ii) - 1.5 + 0.001);
        alpha(ii) = alpha(ii) - tau(ii);
        beta(ii) = beta(ii) - tau(ii);
    else
        tau(ii) = 0;
    end
    if (alpha(ii) > pi)
        alpha(ii) = -2*pi + alpha(ii); 
    end
    if (beta(ii) > pi)
        beta(ii) = -2*pi + beta(ii); 
    end
end

% display block
%plot([1:N], blkRange, 'r--s');
else
    p = 1;
    alpha = -pi+10^-5;
    beta = pi-10^-5;
    tau = 0;   
end
%% Atomic Norm Minimization with block prior
% Other constants
a = tan(alpha/2);
b = tan(beta/2);
d0 = -((a.*b)+1)/2;
d1 = ((1-(a.*b))/4) + (1i*((a+b)/4));

% Complement of J
Jcomp = [1:n]';
Jcomp(J) = [];

% Frequency Recovery by solving dual problem
t = @(ll) ((exp(-1i*tau(ll)*(0:n-1))).');
        
cvx_begin 
variable q(n) complex;
variable Q_1(n, n, p) hermitian;
variable Q_2(n-1, n-1, p) hermitian;
dual variable U{p};

maximize (real(xstar(J)'*q(J)))

subject to
q(Jcomp) == 0;
for ll = 1:p
    [1 (q.*t(ll))'; q.*t(ll) Q_1(:,:,ll)] == semidefinite(n+1,n+1): U{ll};
    Q_2(:,:,ll) == semidefinite(n-1);
    sum(diag(Q_1(:,:,ll),0)) +  sum(d1(ll)*diag(Q_2(:,:,ll),1)) + sum(d0(ll)*diag(Q_2(:,:,ll),0)) + sum(d1(ll)'*diag(Q_2(:,:,ll),-1)) == 1;
    for kk=1:n-1
        sum(diag(Q_1(:,:,ll),-kk)) +  sum(d1(ll)*diag(Q_2(:,:,ll),-kk+1)) + sum(d0(ll)*diag(Q_2(:,:,ll),-kk)) + sum(d1(ll)'*diag(Q_2(:,:,ll),-kk-1)) == 0;
    end
end
cvx_end

%% Van Dec 
% Vandermonde decomposition method introduced in Yang's paper is used for simplicity.
famp=[];
amp=[];
for ii=1:p
    u=U{ii}(2,2:n+1).';
    [fEst1, amp1] = VanDec(u);
    famp=[famp;fEst1,amp1];
end
fSort=-sortrows(-famp,2);
fEst=fSort(1:p,1);

%% display Mean Square Error (MSE)
nfEst = size(fEst,2);
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
fprintf('BANM runtime/iteration: %f\n',toc(tBANM));
fprintf('BANM: Euclidean Distance of freq = %f\n', EDFreq);
fprintf('---------------------------------------------------------\n');

end









