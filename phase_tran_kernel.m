function SuccRateMat = phase_tran_kernel(svec,mvec,nc,solver,TN)
% This file is to compute the phase transition curve
% or phase transition figure using Hankel matrix recovery approach 
% or atomic norm minimization approach with fixed frequency separation.
% 
% Inputs
% - svec: number of atoms, either a scalar or an array. Each element should be greater than 1. 
% - mvec: number of measurements, either a scalar or an array. Each
%   element must be within {1,2,...,2*nc-1}
% - nc: 2*nc-1 is the number of total time domain samples.
% - solver: choose the solver for recovery, either 'HMC' or 'ANM'.
% 
% Return
% - SuccRateMat: a 2D array, each element is the successful recovery rate.
%   The row index is the number of atoms, and the column index is the
%   number of measurements
%
% Variables 
% - delta_f: minimal frequency separation.
% 
% Created by JYI, 10/25/2019.
%
%% Parameter setup
trialnum = TN;
delta_f = 1*10^(-4); tol = 10^(-6);
n = 2*nc-1; t = 0:(n-1);

%% Obtain successful recovery rate
sn = numel(svec); mn = numel(mvec);
SuccRateMat = zeros(sn,mn);

for i=1:sn
    s = svec(i);
    for j=1:mn
        m = mvec(j);
        count = 0;
        trial = 0;
        while trial < trialnum
            
            % Prompt
            fprintf('Atom #: %d, Measurements #: %d, Trial: %d/%d\n',s,m,trial+1,trialnum);
            
            % Generate signal
            f(1:s-1) = rand(1,s-1);
            f(s)=f(1)+delta_f;
            c_mag = (1+10.^(0.5*rand(1,s)))';
            c_pha = 2*pi*rand(1,s)';
            c = c_mag.*exp(sqrt(-1)*c_pha);
            V = exp(1i*2*pi*kron(t',f));
            x = V*c; 
            
            Omega=randsample(n,m);
            
            % Solve optimization
            if strcmp(solver,'ANM')
                cvx_begin sdp quiet
                variable u(n,1) complex
                variable ua(n,1) complex
                variable sv complex
                Q = toeplitz(ua);

                minimize real(trace(Q)/n + trace(s))
                subject to
                     [sv u';u Q] >= 0;
                     u(Omega) == x(Omega);

                cvx_end
                clear cvx_*
            elseif strcmp(solver,'HMC')
                
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
                clear cvx_*
                
            end
            
            % Evaluation
            res_err = norm(u-x)/norm(x);
            if res_err < tol
                count = count +1;
            end
            
            trial = trial+1;
            
        end
        
        SuccRateMat(i,j) = count/trialnum;
        
    end
end

end
