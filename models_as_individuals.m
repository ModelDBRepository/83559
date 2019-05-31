function models_as_individuals(n_batches,n_models,structures,n_cells_per_structure,...  
                                extract_thresh,pathroot,pars_file,flags_file,exp_name,type,varargin)

% MODELS_AS_INDIVIDUALS simulate a batch of models
%   MODELS_AS_INDIVIDUALS(B,M,S,N,T,PATH,PARS,FLAGS,NAME,TYPE)
%   run a batch of models-as-animals simulation protocol, specifying: B
%   batches of M models each, sampling N cells from structures specifed as
%   strings in cell array S (array N must be the length of the structure
%   cell array) - each cell having a mean firing rate of at least T
%   spikes/s. Results and analysis stored in location PATH, using parameter
%   file PARS and flags file FLAGS, saving under experiment name NAME. The
%   analysis files used are defined by TYPE:
%       'SG' - analyse STN-GPe cell properties (uses batch_analyse_stngoe
%       and combine_stngpe_analysis)
%       'sel' - analyse selection and switching properties: will call
%       function to run batch of these and then classify outputs; [may also
%       need to extract spike trains and run analysis to get smoothed individual
%       trains... LATER]
%
%   MODELS_AS_INDIVIDUALS(...,ENG) where ENG is a string specifying the
%   main-file to be called, which handles the network construction and
%   simulation engine. Default is BATCH_BG_heterogenous_AMPA_NMDA
%
%   Mark Humphries 11/05/2006

batch_analysis_list = {};

if nargin >= 11 engine_type = varargin{1}; end

% run whole loop
for loop1 = 1:n_batches
   analysis_list = {}; 
   extract_list = {};
   batch_name = [exp_name num2str(loop1)];  % number each batch within the experiment
   for loop2 = 1:n_models
       fprintf('\n Running model %d of batch %d \n',loop2,loop1);
       % run model 
       r_seed = loop2+(n_models*(loop1-1));        % seed is different for each batch + model combination
       
       if findstr(type,'sel')
            % running selection experiment - pass seed to handling function....   
            [sel_results_list,a_fname,txt_fname] = batch_selection_grid_DA(batch_name,pars_file,r_seed,pathroot,exp_name);
            extract_list{loop2} = sel_results_list;
       else
           if ~exist('engine_type')
               % use default solution engine and main file
               r_fname = BATCH_BG_heterogenous_AMPA_NMDA(1,pars_file,r_seed,pathroot,exp_name);
           else
               eval(['r_fname = ' engine_type '(1,pars_file,r_seed,pathroot,exp_name);']);
           end
           %% extract spike data to use, save to file - datestamp
        
           % extract spikes, saves to unique file 
           % (note that same neuron indices extracted for a given model+batch
           % combination because of use of r_seed)
           e_fname = extract_spikes(r_fname,pathroot,structures,n_cells_per_structure,extract_thresh,r_seed,batch_name);	 
           extract_list{loop2} = e_fname;

           % delete the simulation results file
           delete(r_fname);

           %% pass to analysis functions, save their output to file - datestamp
           ix = findstr(e_fname,'extracted_results');
           a_fname = [e_fname(1:ix-1) 'analysis'];
       end
       analysis_list{loop2} = a_fname;

       % analysis function here
       if findstr(type,'SG')
            batch_analyse_stngpe(e_fname,flags_file,pathroot,a_fname);
       end
   end

  %%% combine analysis results, save results - datestamp
  if findstr(type,'SG')
      c_fname = combine_stngpe_analysis(analysis_list,flags_file,pathroot,batch_name);
  else 
      c_fname = [];
  end
  
  % store batch list
  batch_analysis_list{loop1,1} = c_fname;
  batch_analysis_list{loop1,2} = analysis_list;
  batch_analysis_list{loop1,3} = extract_list;
  
end

% save all filenames resulting from this experiment
time_now = clock;
unique_name = datestr(time_now,30); 

exp_fname = [pathroot exp_name '_' unique_name '_batch.mat'];
save(exp_fname, 'batch_analysis_list')


