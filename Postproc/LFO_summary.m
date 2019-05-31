% script to assess LF0 data-fitting
clear all

%%%%%%%%%%%%%%%% URETHANE COMPLETE MODEL %%%%%%%%%%%%%%%%%%%%%%
%% paths
%batchA_path = '../ResultsArchive/LFO-urethane/CompleteModel/ConditionA/';
%batchB_path = '../ResultsArchive/LFO-urethane/CompleteModel/ConditionB/';
%batchC_path = '../ResultsArchive/LFO-urethane/CompleteModel/ConditionC/';
%batchD_path = '../ResultsArchive/LFO-urethane/CompleteModel/ConditionD/';

% Collaterals, no STN 3rd DA: batch of 50 - Normal model
% THIS SET IS THE CONTROL MODEL FOR THE PAPER
%batchA = 'LFO_a_20060407T004326_batch.mat';
%batchB = 'LFO_b_20060407T011449_batch.mat';
%batchC = 'LFO_c_20060407T003124_batch.mat'; 
%batchD = 'LFO_d_20060407T000419_batch.mat';

%%%%%%%%%%%%%%% URETHANE - NO COLLATERALS %%%%%%%%%%%%%%%%%%%%%%%
% % paths
%batchA_path = '../ResultsArchive/LFO-urethane/NoCollaterals/ConditionA/';
%batchB_path = '../ResultsArchive/LFO-urethane/NoCollaterals/ConditionB/';
%batchC_path = '../ResultsArchive/LFO-urethane/NoCollaterals/ConditionC/';
%batchD_path = '../ResultsArchive/LFO-urethane/NoCollaterals/ConditionD/';
% 
% % No Collaterals, no STN 3rd DA - comparison model FOR PAPER
%batchA = 'LFO_5_4a_20060407T005141_batch.mat';
%batchB = 'LFO_5_4b_20060407T015252_batch.mat';
%batchC = 'LFO_5_4c_20060407T003244_batch.mat';
%batchD = 'LFO_5_4d_20060407T004035_batch.mat';

%%%%%%%%%%%%%%%% Complete Model: NO URETHANE %%%%%%%%%%%%%%%%%%%%%%
% % paths
%batchA_path = '../ResultsArchive/LFO-urethane/NoUrethane/ConditionA/';
%batchB_path = '../ResultsArchive/LFO-urethane/NoUrethane/ConditionB/';
%batchC_path = '../ResultsArchive/LFO-urethane/NoUrethane/ConditionC/';
%batchD_path = '../ResultsArchive/LFO-urethane/NoUrethane/ConditionD/';
 
%batchA = 'LFO_5_5a_20060407T011020_batch.mat';
%batchB = 'LFO_5_5b_20060407T012628_batch.mat';
%batchC = 'LFO_5_5c_20060407T005832_batch.mat';
%batchD = 'LFO_5_5d_20060407T001005_batch.mat';

%%%%%%%%%%%%%%%% URETHANE - CA2+ Mechanism omission %%%%%%%%%%%%%%%%%%%%%
batchA_path = '';
batchB_path = '';
batchC_path = '';
batchD_path = '../ResultsArchive/LFO-urethane/NoCa/ConditionD/';
batchD = 'LFO_5_6d_20060407T001016_batch.mat';

%%%%%%%%%%%%%%%% URETHANE - DA MECHANISM ALTERATIONS %%%%%%%%%%%%%%%%%%%%%%

% % No STN or GP DA mechanism
%batchA_path = '../ResultsArchive/LFO-urethane/NoSTNGP_DA/ConditionA/';
%batchB_path = '../ResultsArchive/LFO-urethane/NoSTNGP_DA/ConditionB/';
%batchC_path = '../ResultsArchive/LFO-urethane/NoSTNGP_DA/ConditionC/';
%batchD_path = '../ResultsArchive/LFO-urethane/NoSTNGP_DA/ConditionD/';
% 
%batchA = 'LFO_5_1a_20060407T004613_batch.mat';
%batchB = 'LFO_5_1b_20060407T011920_batch.mat';
%batchC = 'LFO_5_1c_20060407T002948_batch.mat';
%batchD = 'LFO_5_1d_20060407T000456_batch.mat';

% No STN DA mechanism
%batchA_path = '../ResultsArchive/LFO-urethane/NoSTN_DA/ConditionA/';
%batchB_path = '../ResultsArchive/LFO-urethane/NoSTN_DA/ConditionB/';
%batchC_path = '../ResultsArchive/LFO-urethane/NoSTN_DA/ConditionC/';
%batchD_path = '../ResultsArchive/LFO-urethane/NoSTN_DA/ConditionD/';
% 
%batchA = 'LFO_5_2a_20060407T004637_batch.mat';
%batchB = 'LFO_5_2b_20060407T070147_batch.mat';
%batchC = 'LFO_5_2c_20060407T061658_batch.mat';
%batchD = 'LFO_5_2d_20060407T000453_batch.mat';

% No GP DA mechanism
%batchA_path = '../ResultsArchive/LFO-urethane/NoGP_DA/ConditionA/';
%batchB_path = '../ResultsArchive/LFO-urethane/NoGP_DA/ConditionB/';
%batchC_path = '../ResultsArchive/LFO-urethane/NoGP_DA/ConditionC/';
%batchD_path = '../ResultsArchive/LFO-urethane/NoGP_DA/ConditionD/'; 

%batchA = 'LFO_5_3a_20060407T004251_batch.mat';
%batchB = 'LFO_5_3b_20060407T011649_batch.mat';
%batchC = 'LFO_5_3c_20060407T003058_batch.mat';
%batchD = 'LFO_5_3d_20060407T000413_batch.mat';
 
%%%%%%%%%%%%%%% KETAMINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batchA_path = '../ResultsArchive/LFO-ketamine/CompleteModel/ConditionA/';
% batchB_path = '';
% batchC_path = '';
% batchD_path = '';
% 
% batchA = 'LFO_k_a_20060112T143950_batch.mat';

%% reference data: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
data_Da_means = [8.1; 18; 5.7; 17];     % STN, GP, STN, GP
data_Da_SD = [2.9; 6.3; 4.3; 8.3];  
data_Da_LFO = [100; 0; 0; 0];

data_noDa_means = [18.9; 21.7; 5.8; 18.7];     % STN, GP, STN, GP
data_noDa_SD = [12.2; 8.6; 4.2; 7.3];  
data_noDa_LFO = [100; 89; 20; 15];

%%% for correlations
data_means = [data_Da_means; data_noDa_means];
data_LFO = [data_Da_LFO; data_noDa_LFO];

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

%% load batch lists find the number of batches
n_batches = [];
if exist('batchA_path') & ~isempty(batchA_path)
    load([batchA_path batchA]); 
    batchA_list = batch_analysis_list;
    [n c] = size(batchA_list);
    n_batches = [n_batches n];
end


if exist('batchB_path') & ~isempty(batchB_path)
    load([batchB_path batchB]);
    batchB_list = batch_analysis_list;
    [n c] = size(batchB_list);
    n_batches = [n_batches n];
end

if exist('batchC_path') & ~isempty(batchC_path)
    load([batchC_path batchC]);
    batchC_list = batch_analysis_list;
    [n c] = size(batchC_list);
    n_batches = [n_batches n];
end

if exist('batchD_path') & ~isempty(batchD_path)
    load([batchD_path batchD]);
    batchD_list = batch_analysis_list;
    [n c] = size(batchD_list);
    n_batches = [n_batches n];
end

n_batches = min(n_batches);

model_Da_means = zeros(4,n_batches);
model_Da_SD = zeros(4,n_batches);
model_Da_LFO = zeros(4,n_batches);
model_noDa_means = zeros(4,n_batches);
model_noDa_SD = zeros(4,n_batches);
model_noDa_LFO = zeros(4,n_batches);
num_sig_fits_means = 0;
num_sig_fits_LFO = 0;
batch_fits_means = zeros(n_batches,2);
batch_fits_LFO = zeros(n_batches,2);

% these only needed for other fit-testing options 2+3
%batch_fit_data = zeros(n_batches,2);
%batch_fit_max = zeros(n_batches,2);
%num_sig_fits = zeros(n_batches,2);
%overall_error_per_condition = zeros(n_batches,2);

for loop1 = 1:n_batches
    % get Control data
    %load([batchA_path batchA_list{loop1,2}{1}]); % load first analysis file
    %n_STN_A = length(mean_STN) * length(batchA_list{loop1,2});      % number of STN cells
        
    if exist('batchA_path') & ~isempty(batchA_path)
        load([batchA_path batchA_list{loop1,1}]) 
        model_Da_means(1,loop1) = STN_Hz;
        model_Da_SD(1,loop1) = std_STN;
        model_Da_LFO(1,loop1) = STN_LFO_count ./ n_STN * 100;

        %n_GPe_A = length(mean_GPe) * length(batchA_list{loop1,2});      % number of STN cells
        model_Da_means(2,loop1) = GPe_Hz;
        model_Da_SD(2,loop1) = std_GPe;
        model_Da_LFO(2,loop1) = GPe_LFO_count ./ n_GPe * 100;
    end
    
    % get DA, no Ctx data
    %load([batchB_path batchB_list{loop1,2}{1}]); % load first analysis file
    %n_STN_B = length(mean_STN) * length(batchB_list{loop1,2});      % number of STN cells
    if exist('batchB_path') & ~isempty(batchB_path)
        load([batchB_path batchB_list{loop1,1}]) 

        model_Da_means(3,loop1) = STN_Hz;
        model_Da_SD(3,loop1) = std_STN;
        model_Da_LFO(3,loop1) = STN_LFO_count ./ n_STN * 100;

        %n_GPe_B = length(mean_GPe) * length(batchB_list{loop1,2}); 
        model_Da_means(4,loop1) = GPe_Hz;
        model_Da_SD(4,loop1) = std_GPe;
        model_Da_LFO(4,loop1) = GPe_LFO_count./ n_GPe * 100;;
    end
    
    % get no DA, Ctx data
    if exist('batchC_path') & ~isempty(batchC_path)
        load([batchC_path batchC_list{loop1,1}]) 
        %load([batchC_path batchC_list{loop1,2}{1}]); % load first analysis file
        %n_STN_C = length(mean_STN) * length(batchC_list{loop1,2});      % number of STN cells

        model_noDa_means(1,loop1) = STN_Hz;
        model_noDa_SD(1,loop1) = std_STN;
        model_noDa_LFO(1,loop1) = STN_LFO_count./ n_STN * 100;

        %n_GPe_C = length(mean_GPe) * length(batchC_list{loop1,2}); 
        model_noDa_means(2,loop1) = GPe_Hz;
        model_noDa_SD(2,loop1) = std_GPe;
        model_noDa_LFO(2,loop1) = GPe_LFO_count./ n_GPe * 100;
    end
    
    % get no DA, no Ctx data
    if exist('batchD_path') & ~isempty(batchD_path)
        load([batchD_path batchD_list{loop1,1}]) 
        %load([batchD_path batchD_list{loop1,2}{1}]); % load first analysis file
        %n_STN_D = length(mean_STN) * length(batchD_list{loop1,2});      % number of STN cells

        model_noDa_means(3,loop1) = STN_Hz;
        model_noDa_SD(3,loop1) = std_STN;
        model_noDa_LFO(3,loop1) = STN_LFO_count./ n_STN * 100;

        %n_GPe_D = length(mean_GPe) * length(batchD_list{loop1,2}); 
        model_noDa_means(4,loop1) = GPe_Hz;
        model_noDa_SD(4,loop1) = std_GPe;
        model_noDa_LFO(4,loop1) = GPe_LFO_count./ n_GPe * 100;
    end
    
    %% compute overall fit to mean firing rate and LFO frequency data!
    
    % create vectors of fitted data
    model_means = [model_Da_means(:,loop1); model_noDa_means(:,loop1)];
    model_LFO = [model_Da_LFO(:,loop1); model_noDa_LFO(:,loop1)];
    
    %%%% option 1: compute separate correlations, using Spearman's rank and
    %%%% doing one-tailed test for direction
    [r_mean,p_mean] = corr(data_means,model_means,'tail','gt','type','Spearman');
    [r_LFO,p_LFO] = corr(data_LFO,model_LFO,'tail','gt','type','Spearman');
    
    batch_fit_means(loop1,:) = [r_mean^2,p_mean];
    batch_fit_LFO(loop1,:) = [r_LFO^2,p_LFO];
    
    %   find significance
    if p_mean < 0.05
        num_sig_fits_means = num_sig_fits_means +1;
    end
    
     if p_LFO < 0.05
        num_sig_fits_LFO = num_sig_fits_LFO +1;
    end
    

    
%     %%%% option 2: rescale to study data
%     % 0.  re-scale study data (could put this outside loop, but left here
%     % for comprehension)
%     data_scaling = max(data_LFO) / max(data_means);
%     data_scaled_means = data_means .* data_scaling;
%     all_data = [data_scaled_means; data_LFO]; 
%     
%     % 1. re-scale model data too, and combine
%     model_means_scaled = model_means .* data_scaling;
%     all_model_data = [model_means_scaled; model_LFO];
%     
%     % 2. do correlation    
%     [r,pval] = corr(all_data,all_model_data,'tail','gt','type','Spearman');
%     
%     batch_fit_data(loop1,:) = [r, pval];
%     
%     % 3. find significance
%     if pval < 0.05
%         num_sig_fits(loop1,1) = num_sig_fits(loop1,1)+1;
%     end
%     
%     % 4. find average error per condition
%     overall_error_per_condition(loop1,1) = sum(abs(all_data - all_model_data)) / length(all_data);
%     
%     %%%% option 3: rescale to maximum mean firing rate in model or data to
%     %%%% ensure that error bounded between 0-100%
%     
%     % 1. re-scale all mean data
%     scaling = max(max(data_LFO),max(model_LFO)) / max(max(data_means),max(model_means));
%     data_scaled_means = data_means .* scaling;
%     all_data = [data_scaled_means; data_LFO]; 
%     model_means_scaled = model_means .* scaling;
%     all_model_data = [model_means_scaled; model_LFO];
%     
%     % 2. do correlation
%     [r_max,p_max] = corr(all_data,all_model_data,'tail','gt','type','Spearman');
%     batch_fit_max(loop1,:) = [r_max, p_max];
%     
%     % 3. find significance
%     if p_max < 0.05
%         num_sig_fits(loop1,2) = num_sig_fits(loop1,2)+1;
%     end
%     
%     % 4. find average error per condition
%     overall_error_per_condition(loop1,2) = sum(abs(all_data - all_model_data)) / length(all_data);
% 
    
    %% save data
%     temp = model_Da_means(:,loop1);
%     save(['batch_' num2str(loop1) '_model_Da_means.txt'],'temp','-ascii')
%     temp = model_Da_SD(:,loop1);
%     save(['batch_' num2str(loop1) '_model_Da_SD.txt'],'temp','-ascii')
%     temp = model_Da_LFO(:,loop1);
%     save(['batch_' num2str(loop1) '_model_Da_LFO.txt'],'temp','-ascii')
%     temp = model_noDa_means(:,loop1);
%     save(['batch_' num2str(loop1) '_model_noDa_means.txt'],'temp','-ascii')
%     temp = model_noDa_SD(:,loop1);
%     save(['batch_' num2str(loop1) '_model_noDa_SD.txt'],'temp','-ascii')
%     temp = model_noDa_LFO(:,loop1);
%     save(['batch_' num2str(loop1) '_model_noDa_LFO.txt'],'temp','-ascii')
end

% save overall matrices of summaries too - useful for SigmaPlot
model_Da_means = model_Da_means';       % transpose so that condition is down column
model_noDa_means = model_noDa_means'; 
model_Da_SD = model_Da_SD'; 
model_noDa_SD = model_noDa_SD'; 
model_Da_LFO = model_Da_LFO';       % transpose so that condition is down column
model_noDa_LFO = model_noDa_LFO'; 

save all_model_DA_means.txt model_Da_means -ascii
save all_model_noDA_means.txt model_noDa_means -ascii
save all_model_DA_SD.txt model_Da_SD -ascii
save all_model_noDA_SD.txt model_noDa_SD -ascii
save all_model_DA_LFO.txt model_Da_LFO -ascii
save all_model_noDA_LFO.txt model_noDa_LFO -ascii

% save fit data
save('batch_fit_means.txt','batch_fit_means','-ascii');
save('batch_fit_LFO.txt','batch_fit_LFO','-ascii');
save('batch_num_sig_fits_means.txt','num_sig_fits_means','-ascii');
save('batch_num_sig_fits_LFO.txt','num_sig_fits_LFO','-ascii');
%% look at problem data - Ctx, no DA GPe LFO count

total_fits = num_sig_fits_means + num_sig_fits_LFO
num_sig_fits_means
num_sig_fits_LFO




