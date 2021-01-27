function [ EDFreq,blkFin ] = BL1M( xstar, J, A, N, k, n, fplot, freq, ck, MaxItr )
%% Block Iterative Reweighted L1 Minimization (BL1M)
%
% - Input
% xstar: all time data samples
% J: randomly chosen time indices
% A: partial DFT matrix 
% N: number of discretization (= number of columns in the matrix A)
% k: number of sparsity
% n: number of signal dimension (n x 1)
% fplot: frequency plot for display purpose
% freq: true frequency in [0,1] (for error calculation)
% ck: true frequency coefficient (display purpose)
% MaxItr: maximum number of iteration
%
% - Output
% EDFreq: Euclean distance between true freq. and the estimated freq.
% blkFin: final block information for later use in BANM-Mix
%
% Michael Myung Cho
% myung-cho@uiowa.edu
%%

%% Initial Setting for block size, weights, Indices
fEstOld = zeros(N,1);
xEstOld = zeros(n,1);
BW = 320;       % Block Width
gridCol=[];     % for adaptive grid
w = ones(N,1);  % weights
Jc = [1:n]';    % Complement of J, (J: chosen time index)
Jc(J) = []; 

%% block predefine for adaptive grid (adpative block size) 
s = [1:N]';
b=zeros(BW+1,N);
for ii=1:N
    t = circshift(s,BW/2-ii+1);
    iS=mod(ii-BW/2,N);
    iE=mod(ii+BW/2,N);
    if iS < iE
        t=[iS:iE]';
    else
        t=[[iS:N]';[1:iE]'];
    end
    b(:,ii)=t;
end
b(find(b==0))=N;
mb=BW/2+1;

%% Iteration 
for itr=1:MaxItr
    tB1LM=tic;
    %% small value for weight to prevent infinity.
    if itr <= 5
        epsilon = 1/2^itr;
    else
        epsilon = 1/2^8;
    end
    
    %% BL1M (primal problem) implemented in CVX 
    % adaptive grid
    if itr==1
        gridCol=[1:2^4:N]'; 
    else
        wMid = (min(w) + max(w))/2;
        Fi = find(w < wMid);
        gridCol=unique([gridCol;Fi]);
    end
    pN=size(gridCol,1); 
    
    % CVX starts
    cvx_quiet true
    cvx_precision default
    
    cvx_begin
        variable ctemp(pN) complex;
        minimize norm(w(gridCol).*ctemp,1);
    subject to
        xstar(J) == A(J,gridCol)*(ctemp);
    cvx_end
    c = zeros(N,1); % frequency coefficient which we want to find by solving optimization problem
    c(gridCol)=ctemp;

    %% calculation
    fEst = c;
    xEst = A*c;
    %% display signal result and histogram result
%     figure(1);
%     plot([1:N], abs(fplot),'b--o', [1:N],abs(fEst),'r--s')
%     figure(3);
%     plot([1:N], w(1:N),'r--s')
    
    %% Stop creteria
    fED=norm(fEst - fEstOld,2);
    if (  fED < 5*10^-4 && itr>=2) || (itr == MaxItr)
        
        % Center mass calculation
        avgF = [];
        wMid = (min(w) + max(w))/2;
        Fi = find(w < wMid);
        tBlk = zeros(N,1);
        tBlk(Fi) = 1;
        for ii=1:N
            sFlag = 0;
            eFlag = 0;
            cir=0;
            if ii==1
                if (tBlk(N) == 1) && (tBlk(ii) == 1)
                    for jj=N:-1:1
                        if(tBlk(jj) == 1 && tBlk(jj-1) == 0)
                            sPoint=jj;
                            sFlag=1;
                            cir=1;
                            break;
                        end
                    end
                    for jj=1:N
                        if(tBlk(jj) == 1 && tBlk(jj+1) == 0)
                            ePoint = jj;
                            eFlag = 1;
                            break;
                        end
                    end
                end
                if (tBlk(N) == 0) && (tBlk(ii) == 1)
                    sPoint=1;
                    sFlag=1;
                    for jj=1:N
                        if(tBlk(jj) == 1 && tBlk(jj+1) == 0)
                            ePoint = jj;
                            eFlag = 1;
                            break;
                        end
                    end
                end
            else
                if (tBlk(ii-1) == 0) && (tBlk(ii) == 1)
                    sPoint = ii;
                    sFlag = 1;
                    for jj=sPoint:N-1
                        if (tBlk(jj)== 1) && (tBlk(jj+1) == 0)
                            ePoint = jj;
                            eFlag = 1;
                            break;
                        end
                    end
                end
            end
            if sFlag == 1 && eFlag ==1
                sumFM = 0;
                if cir == 0
                    EstFM = sum(abs(fEst(sPoint:ePoint)));
                    for jj=sPoint:ePoint
                        sumFM = sumFM + jj*abs(fEst(jj,1));
                    end
                else
                    EstFM = sum(abs(fEst(sPoint:N)))+ sum(abs(fEst(1:ePoint)));
                    for jj=1:ePoint
                        sumFM = sumFM + jj*abs(fEst(jj,1));
                    end
                    for jj=sPoint:N
                        sumFM = sumFM + jj*abs(fEst(jj,1));                        
                    end
                end

                avgFi = floor(sumFM/EstFM);
                if avgFi > N
                    avgFi = N - avgFi;
                else
                    avgF = [avgF;avgFi];
                end
            end
        end
        
        % display Mean Square Error (MSE)
        nfEst = size(avgF,1);
        fErr = 0;
        if nfEst >= k
            freqC = [freq;min(freq)+1;max(freq)-1];
            for ii=1:nfEst
                % freq error in circle
                [res, ~] = sort(abs(freqC - avgF(ii)/N),'ascend');
                fErr = fErr + res(1)^2;
            end
        else % case: nEstF < k
            avgFC = [avgF;min(freq)+N;max(freq)-N];
            for ii=1:k
                % freq error in circle
                [res, ~] = sort(abs(freq(ii) - avgFC./N),'ascend');
                fErr = fErr + res(1)^2;
            end
        end
        EDFreq = sqrt(fErr);
        fprintf('B1LM runtime/iteration: %f (dis:%f)\n', toc(tB1LM),fED);
        fprintf('BL1M_primal: Euclidean Distance of freq = %f (itr = %d)\n', EDFreq,itr);
        fprintf('---------------------------------------------------------\n');
        blkFin = zeros(1,N);
        bInd = find(w < wMid);
        blkFin(bInd) = 1;
        break;
    end
    fEstOld = fEst;
    xEstOld = xEst;

    %% block circular
    % calculating weights
    for ii=1:N
        w(ii) = 1/(sum(abs(c( b(mb-BW/2:mb+BW/2,ii) ))) + epsilon);
    end  
    % adaptive grid (reducing the block width)
    if itr < 5
        BW=BW/2;
    end
fprintf('B1LM runtime/iteration: %f (dis:%f)\n', toc(tB1LM),fED);
end




end

