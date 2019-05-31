% script to combine categories from selection and switching experiments
% Mark Humphries 31/1/2006

clear all

% normal model
%batch_name = 'selection_20060407T015909_batch.mat';
%batch_path = '../ResultsArchive/Selection/Normal/';

% low DA model
%batch_name = 'selection_lowDA_20060407T020953_batch.mat';
%batch_path = '../ResultsArchive/Selection/LowDA/';

% high DA model (da = 0.8?)
batch_name = 'selection_highDA_20060407T025808_batch.mat';
batch_path = '../ResultsArchive/Selection/HighDA/';

%% set paths for interactive sessions
[a host] = system('hostname');

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

% get batch analysis lists
load([batch_path batch_name]);
[r n_models] = size(batch_analysis_list{2});

distribution = zeros(n_models,5);
% loop over all models, load analysis 
for loop = 1:n_models
   load([batch_path batch_analysis_list{2}{loop}]) 
   categories = batch_sum(:,5:end);
   distribution(loop,:) = sum(categories) ./ 4;      % where 4 is the bubble scale
end

mean_categories = mean(distribution);
std_categories = std(distribution);

save selection_batch_categories.txt distribution -ascii
save selection_mean_categories.txt mean_categories -ascii
save selection_std_categories.txt std_categories -ascii




