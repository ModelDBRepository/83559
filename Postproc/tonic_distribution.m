%%% script to look at tonic rate distributions
% Mark Humphries 2/2/2006

clear all

batch_path = '../ResultsArchive/tonic/WithoutCollaterals/'; 
batch_name = 'tonic_20060407T133735_batch.mat';   % without collaterals

%batch_path = '../ResultsArchive/tonic/';
%batch_name = 'tonic_20060406T192735_batch.mat';     % with collaterals in GP/SNr

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

n_batches = length(batch_analysis_list);

STN_array = zeros(n_batches,2);
GPe_array = zeros(n_batches,2);
GPi_array = zeros(n_batches,2);

for loop = 1:n_batches
    load([batch_path batch_analysis_list{loop,1}]);
    STN_array(loop,1) = STN_Hz;
    STN_array(loop,2) = sem_STN;
    GPe_array(loop,1) = GPe_Hz;
    GPe_array(loop,2) = sem_GPe;
    GPi_array(loop,1) = GPi_Hz;
    GPi_array(loop,2) = sem_GPi;
end

batch_mean_STN = mean(STN_array(:,1));
batch_mean_GPe = mean(GPe_array(:,1));
batch_mean_GPi = mean(GPi_array(:,1));
batch_se_mean_STN = std(STN_array(:,1)) ./ sqrt(n_batches);
batch_se_mean_GPe = std(GPe_array(:,1)) ./ sqrt(n_batches);
batch_se_mean_GPi = std(GPi_array(:,1)) ./ sqrt(n_batches);

cl_fig
figure
hist(STN_array(:,1),10);
title('Distribution of STN means')

figure
hist(GPe_array(:,1),10);
title('Distribution of GP means')

figure
hist(GPi_array(:,1),10);
title('Distribution of SNr means')

tile

%% save data
save STN_tonic_batch.txt STN_array -ascii
save GPe_tonic_batch.txt GPe_array -ascii
save GPi_tonic_batch.txt GPi_array -ascii




