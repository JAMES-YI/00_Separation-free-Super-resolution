%% Frequency recovery from disretization 
% - Abstract: 
% finding frequency information with a few randomly chosen time sampled data
% 
% - Paper: 
% Block Iterative Reweighted Algorithms for Super-Resolution of Spectrally Sparse Signals
% IEEE Signal Processing Letters, vol. 22, no. 12, pp. 2319–2323, 2015.
% 
% - Abbreviation
% Previous Algorithms
% 1. SANM: Standard Atomic Norm Minimization 
%   (Paper: Compressed Sensing off the Grid by G. Tang, et al.) 
% 2. RANM: Reweighted Atomic Norm Minimization 
%   (Paper: Enhancing Sparsity and Resolution via Reweighted Atomic Norm Minimization by Z. Yang, et al.)
% Proposed Algorithms
% 1. BANMx: Block Iterative Reweighted Atomic Norm - Mix. 
% 2. BL1M: Block Iterative Reweighted L1 Minimization
%
% Myung (Michael) Cho 
% University of Iowa 
% Sept. 11, 2014
% Contect: myung-cho@uiowa.edu
%%
clear; 
clc;

%% Dimensions of the problem
n = 32;           % Number of time-sample (=number of rows in sensing matrix A)
N = 2^14;         % Number of columns in sensing matrix (discretization) 2^13=8192, 2^14=16384
ks = 5;           % Number of non-zero elements (number of frequencies)
randSamples = 20; % random time sample

%% Generate dicretized sensing matrix
A = generationFourierMatrixAtomByAtom( n, N );
 
%% iteration for testing start
% Initialize the number of recovered frequencies for each iteration
res =[];

%% save data to xls file
countXLS=1;    
clk=clock;
datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5))]);
filename = strcat('(',num2str(n),'x',num2str(N),')','Prob-SANM-RANM-BANMMix-BL1M',datestmp,'.xlsx');
col_name = {'n','N','k','randSample','nSuccSANM','nSuccRANM','nSuccBANMx','nSuccBL1M';};
xlswrite(filename,col_name);

%% iteration
nTrial = 50; % number of trials
MaxItr = 20; % Maximum iteration number
for k = ks
    for randSample = randSamples
        randSample
        nSuccBL1M = 0; 
        nSuccSANM = 0;
        nSuccBANMx = 0;
        nSuccRANM = 0;
        for iTrial=1:nTrial
            iTrial
            %% Create a signal with s different frequencies 
            scaling_flag = 0; % scaling_flag = 1 will scale-down/divide the signal by sqrt(|J|)
            noise_flag = 0;   % no noise is added
            WL = 0;
            WH = 1;
            % off-the-grid orignal signal case
            [xstar, freq, phases, ck] = generateBandLimitedSignalOffGrid(n, k, noise_flag, scaling_flag, WL, WH, 'random_anysep');
            % Find maximum frequency separation
            freq_sort = sort(freq', 'ascend');
            freq_diff = [diff(freq_sort) (freq_sort(1) + (1-freq_sort(end)))];
            delta_f = min(freq_diff);
            % display of original signal
            gridfreq=floor(freq*N);
            fplot=zeros(1,N);
            fplotFreq=zeros(1,N);
            fplot(1,gridfreq + N*(gridfreq == 0))=ck.';
            
            %% Select samples uniformly at random in the interval
            unif_rand_samples = randperm(n);
            T = (1:n)';         % total time sample
            J = T(unif_rand_samples(1:randSample));
            J = sort(J);        % randomly chosen time sample

            %% BL1M
            tBL1M = tic;
            [ED_F_BIT, blkFin] = BL1M( xstar, J, A, N, k, n, fplot, freq, ck, MaxItr );
            tBL1M2 = toc(tBL1M);
            if (ED_F_BIT < 10^-3)
                 nSuccBL1M = nSuccBL1M + 1;
            end
            %% RANM            
            % In order to run SANM, RANM code is used by setting the
            % MaxItr number to 1.
            tSANM = tic;
            [ED_F_SANM ]  = RANM_dual(xstar, J, n, k, fplot, freq, ck, 1); 
            tSANM2 = toc(tSANM); 
            if (ED_F_SANM < 10^-3) 
                 nSuccSANM = nSuccSANM + 1;
            end            
             
            %% BANM-Mix.
            tBANMn = tic;
            [ ED_F_BANMx ] = BANM_Mix( xstar, J, N, k, n, fplot, freq, ck, blkFin, delta_f );
            tBANM2n = toc(tBANMn);
            if (ED_F_BANMx < 10^-3 || ED_F_BIT < 10^-3)
                 nSuccBANMx = nSuccBANMx + 1;
            end
            
            %% RANM
            tRANM = tic;
            [ ED_F_RANM ] = RANM_dual(xstar, J, n, k, fplot, freq, ck, MaxItr);
            tRANM2 = toc(tRANM);
            if (ED_F_RANM < 10^-3)
                 nSuccRANM = nSuccRANM + 1;
            end
        end
        resBuf = [n,N,k,randSample, nSuccBL1M, nSuccBANMx, nSuccSANM, nSuccRANM ];
        %% save data to xls file
        countXLS = countXLS+1;
        xlRange = sprintf('A%d',countXLS+1);
        xlsInsert = resBuf;
        xlswrite(filename, xlsInsert, 'Sheet1', xlRange);
    end
end










