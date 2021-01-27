function [xstar, freq, phases, ck] = generateBandLimitedSignalOffGrid(n, s, noise_flag, varargin)
%  Kumar Vijay Mishra
%  The University of Iowa, Iowa City, IA
%  Generate bandlimited frequency-sparse regular signal vector
%  Sep 24 2013
%
%  Inputs:
%  n = number of samples
%  s = number of frequencies
%  noise_flag = Add noise if this is equal to 1
%  varargin are 
%   scaling_flag = If this flag is 1, xstar will be divided by sqrt(|J|)
%   WL = lower limit of the bandwidth (default 0)
%   WH = lower limit of the bandwidth (default 1)
%   spacing = 'equi_spaced' - s frequencies picked from equi-divided [0 1]
%             'random_minsep' - s frequencies randomly picked from [0 1]
%                               with minimum specified separation
%             'random_anysep' - s frequencies randomly picked from [0 1]
%   delta_f = frequency spacing (used only for random_minsep option)
%
%  Outputs:
%  xstar = Regular signal vector
%  freq = True frequencies 
%  ck = True amplitudes
%  phases = True phases

time = (0:(n-1))'; % Time
w =  randn(s, 1); % zero mean unit variance Gaussian random variable
ck = sqrt(0.5 + (w.^2)); % Fading amplitudes
%ck = ck./max(ck);

if (isempty(varargin))
    scaling_flag = 0;
    WL = 0;
    WH = 1;
    spacing = 'random_minsep';
    delta_f = 1/((n-1)/4);
elseif (length(varargin) == 1)
    scaling_flag = varargin{1};
    WL = 0;
    WH = 1;
    spacing = 'random_minsep';
    delta_f = 1/((n-1)/4);
elseif (length(varargin) == 2)
    scaling_flag = varargin{1};
    WL = varargin{2};
    WH = 1;
    spacing = 'random_minsep';
    delta_f = 1/((n-1)/4);
elseif (length(varargin) == 3)
    scaling_flag = varargin{1};
    WL = varargin{2};
    WH = varargin{3};
    spacing = 'random_minsep';
    delta_f = 1/((n-1)/4);
elseif (length(varargin) == 4)
    scaling_flag = varargin{1};
    WL = varargin{2};
    WH = varargin{3};
    spacing = varargin{4};
    delta_f = 1/((n-1)/4);
elseif (length(varargin) >= 5)
    scaling_flag = varargin{1};
    WL = varargin{2};
    WH = varargin{3};
    spacing = varargin{4};
    delta_f = varargin{5};    
end

% Pick s freqs with the specified spacing type in the interval [WL WH]
freq = pickFrequencies(s, delta_f, WL, WH, spacing);

% Get phases - they can either be Bernoulli or Uniform
phases = 2*pi*rand(s, 1); % Uniform distribution on unit circle

% Initialize the frequency sparse signal
xstar = zeros(n, 1);

% Generate the frequency sparse signal 
 %freq = [0.3112    0.2542     0.1895]
for k = 1:s
    xstar = xstar + abs(ck(k))*exp(1j*((2*pi*freq(k)*time) + phases(k))); 
end

% Add complex noise
if noise_flag == 1
    noise_bound = 2;
    noise_ang = rand(1,n)*(2*pi);
    noise_mag = sqrt(rand(1,n))*noise_bound;
    n_complex = (noise_mag.*cos(noise_ang)) + 1i*(noise_mag.*sin(noise_ang));
    xstar = xstar + n_complex';
end

% Scaling: Divide the signal model by the sqrt of the cardinality of the set. Here
% |J| = n
if scaling_flag == 1
    xstar = xstar/(sqrt(n));
end
%************************************************************************

