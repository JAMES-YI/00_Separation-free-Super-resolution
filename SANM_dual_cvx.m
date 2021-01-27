function [rx,err]=SANM_dual_cvx(ox,K,N,kk,stepsize,epsilon)
%% Standard Atomic Norm Minimization
KC = (1:N)';
KC(K) = [];

% solve the dual problem
cvx_quiet true
cvx_precision default
cvx_solver sdpt3

cvx_begin sdp

  variable X(N+1,N+1) hermitian;
  dual variable U;
  
  X >= 0 : U;
  X(1+KC,1) == 0;
  for j = 0:N-1
    sum(diag(eye(N) - X(2:N+1,2:N+1), j)) == 0;
  end
  X(1,1) == 1;

  minimize real( trace( ox(K)'*X(K+1,1)))
  
cvx_end
% Y estimate
rx = U(2:1+N,1)*2;

err = norm(rx-ox)/norm(ox);
end