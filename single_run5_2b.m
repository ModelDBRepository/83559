function single_run5_2b()

% run a single model and postprocess it

[a host] = system('hostname');
fprintf(['\n' host '\n']);

%%% set requisite paths!
if findstr(host, 'iceberg'); % on iceberg
     fprintf('\n On ICEBERG \n');
     system_os = 'unix';
     ice_path1 = genpath('/home1/pc/pc1mdh/BG spiking model');
     ice_path2 = genpath('/home1/pc/pc1mdh/Matlab Tools');
     path(path, ice_path1);
     path(path, ice_path2);
elseif (findstr(host, 'node') | findstr(host,'ace')) % on ACE
     system_os = 'unix';
     ace_path1 = genpath('/home/mark/SpikingModel');
     path(path, ace_path1);
     fprintf('\n On ACE \n');
else
     system_os = 'xp';
     fprintf('\n On XP \n');
end

%%%% always need all of these parameters %%%%%%
n_batches = 50;
n_models = 7; % 7
structures = {'STN' 'GPe'};
n_cells_per_structure = [2 4];
extract_thresh = 0;

pathroot = ['/home/mark/SpikingModel/ResultsArchive/LFO-urethane/NoSTN_DA/ConditionB/'];  % path to save files
pars_file = 'pars5_2b';
flags_file = 'sum_flags45';

exp_name = 'LFO_5_2b';
type = 'SG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
models_as_individuals(n_batches,n_models,structures,n_cells_per_structure,...  
                                extract_thresh,pathroot,pars_file,flags_file,exp_name,type)
                            
toc
