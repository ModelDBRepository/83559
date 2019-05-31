function e_fname = extract_spikes(r_fname,pathroot,structures,n_cells,threshold,seed,varargin)

% EXTRACT_SPIKES samples spike trains from results
%   EXTRACT_SPIKES(R,PATH,S,N,T,SEED) extracts N_i spikes trains from the ith nucleus named in cell array S
%   (where the names must match the indices arrays e.g. SD1, STN etc) from the results file supplied
%   in string R, on the path PATH. The spike trains selected 
%   will have a mean firing rate (over entire period) of at least T. (N is thus a vector). The 
%   random number generator is initiated with SEED - this allows for
%   recallable random cell selection!
%
%
%   EXTRACT_SPIKES(..,P,FLAG) where P is a string creating a meaningful
%   prefix for the extracted results file
%
%   Set options using FLAG: 
%       'c' -  NOT IMPLEMENTED YET adding this flag ensures that N spikes trains are sampled from each 
%             channel in each structure, and are saved in accordingly named arrays (set to '' to omit). 
%  
%
%   Returns the name of the file produced.
%
%   Mark Humphries 4/1/2006

prefix = '';

rand('state',seed);

if nargin >= 7 
    if ~isstr(varargin{1})
        error('Supplied filename prefix is not a string')
    else
        prefix = varargin{1};
    end
end

% load the simulation results file
load(r_fname)

t_out_t = double(out_t) .* dt; % convert to time

n_structures = length(structures);

if n_structures ~= length(n_cells)
    error('Cells-per-structure and structure arrays do not match');
end


time_now = clock;
unique_name = datestr(time_now,30);  

e_fname = [prefix '_' unique_name '_extracted_results'];
path_fname = [pathroot e_fname '.mat'];

for loop = 1:n_structures
    idxs = structures{loop};
     
   
    for loop2 = 1:n_cells(loop)
        has_spikes = 0; 
        eval(['n_possibles = length(' idxs ');']); 
        rand_idx = randperm(n_possibles);
        counter = 1;    % counter for the permuted sequence - if cell already discarded, no need to test it again!
        
        while ~has_spikes 
           eval(['mean_rate = sum(out_n==' idxs '(rand_idx(counter))) / time_seconds;']);
           counter = counter + 1;  
           if mean_rate >= threshold has_spikes = 1; end
           if counter > n_possibles error('Not enough cells had firing rates that exceeded the extraction threshold'); end            
        end
        eval(['ts = t_out_t(out_n == ' idxs '(rand_idx(counter)));']);
        eval([idxs '_times{loop2} = ts;']);     % generate cell array of spike trains for use in the analysis routines
    end
    
    % save newly created cell arrays
    if loop > 1
        eval(['save(path_fname,''' idxs '_times'',''-append'');'])
    else
        eval(['save(path_fname,''' idxs '_times'');']);
    end
end




% also save the necessary simulation data (including the input arrays)
save(path_fname,'input_array','trace_vals','trace_n',...
        'R','theta','tau_m','tau_AMPA','tau_NMDA','tau_GABAa','spon','mlimit',...
        'dt','time_seconds','p_connect',...
        'n_channels','n_neurons','n_sources','n_inputs','neurons_per_nucleus','neurons_per_channel', 'n_nuclei',...
        'delays','proportions','link_n','link_wAMPA','link_wNMDA','link_wGABAa','link_t',...
        't_on', 't_off', 'step_size', 'switches',...
        'GPi','STN','GPe','EXT','SD1','SD2','t_out_t','out_n','in_n','in_t','-append');
    
   if exist('fast_weights')
       save(path_fname,'fast_weights','slow_weights','-append');
   else
       save(path_fname,'weights','-append');
   end
        
