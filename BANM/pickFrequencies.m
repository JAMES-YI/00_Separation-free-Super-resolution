function freq = pickFrequencies(s, delta_f, WL, WH, spacing)
%  Kumar Vijay Mishra
%  The University of Iowa, Iowa City, IA
%  Generate input regular signal vector for off-the-grid algorithm
%  Sept 23 2013
%
%  Inputs:  
%  s = number of frequencies
%  delta_f = separation between frequencies
%  WL = lower limit of the bandwidth (default 0)
%  WH = lower limit of the bandwidth (default 1)
%  spacing = 'equi_spaced' - s frequencies picked from equi-divided [0 1]
%             'random_minsep' - s frequencies randomly picked from [0 1]
%                               with minimum specified separation
%             'random_anysep' - s frequencies randomly picked from [0 1]
%
%  Outputs:
%  freq = True frequencies 

if strcmp(spacing, 'equi_spaced')
    
    % Generate s equispaced frequencies in [0 1]
    % First pick a random start_freq
    start_freq = WL + (WH-WL).*rand(1,1);
    % Then pick other s-1 equi-spaced frequencies
    translated_freqs = linspace(start_freq, start_freq+(WH-WL), s+1);
    freq = mod(translated_freqs, WH-WL);    
    freq = WL + freq(1:s); % Eliminate the repeated (s+1)th frequency
    
elseif strcmp(spacing, 'random_minsep')
    
    % Get s frequencies randomly selected from a uniform distribution U~[0 1]
    % with a constraint that their separation should be delta_f    
    % UPDATE [Sep 20, 2013]: Following can generate even more random 
    % frequencies with delta_f separation in the range [0 1]
    % Try a maximum of 20 times if there is an error
    max_trials = 20;
    num_trials = 0;
    freq_not_found = true;
    while(freq_not_found)
        freq = random_minsep(s, delta_f, WL, WH);
        if (num_trials <= max_trials)
            if isempty(freq) == 1
                num_trials = num_trials + 1;
            else
                freq_not_found = false;
            end
        else
            error(['Error: Failed to generate frequencies after ' num2str(max_trials) ...
                   ' attempts. Either s or delta_f are bigger. Consider revising the values.']);        
        end
    end % end of while(freq_not_found) loop
    
else
    
    % Generate frequencies randomly with any separation
    freq = WL + (WH-WL).*rand(s,1);
    
end

% Sort all the frequencies before returning to the main program
freq = sort(freq, 'descend');

