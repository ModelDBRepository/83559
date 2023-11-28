% script to get Gamma-osc results for paper, based on models-as-animals
% simulation results
%
% Mark Humphries 31/5/2006

clear all

control_batch = 'Gamma_a_20060407T204507_batch.mat';
%d2_batch = 'Gamma_b_20060407T152129_batch.mat';
% treat NMDA manipulation as equivalent to D2 agonist
d2_batch = 'Gamma_noNMDAinGP_a_20060512T211003_batch.mat'; 

control_path = '../ResultsArchive/Gamma-band/CompleteModel/ConditionA/';
%d2_path = '../ResultsArchive/Gamma-band/CompleteModel/ConditionB/';
% treat NMDA manipulation as equivalent to D2 agonist
d2_path = '../ResultsArchive/Gamma-band/NoNMDAinGP/';

analyse = 'STN'; % either 'STN' or 'GP'

x = 0:5:100;            % bin edges for peaks and firing

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

%% load in control and D2-case lists
load([control_path control_batch]);
control_list = batch_analysis_list;

load([d2_path d2_batch]);
d2_list = batch_analysis_list;

[n_control_batches c] = size(control_list);
[n_d2_batches c] = size(d2_list);

%% load first batch control combined analyses to get no. units etc 
n_control_models = length(control_list{1,2});
n_d2_models = length(d2_list{1,2});

switch analyse
    case 'STN'
        load([control_path control_list{1,2}{1}],'mean_STN','STN_ISIhist')    
        ISI_hist_length = length(STN_ISIhist{1});
        n_control = length(mean_STN) * n_control_models;       % total number of cells sampled
        
        load([d2_path d2_list{1,2}{1}],'mean_STN')
        n_d2 = length(mean_STN) * n_d2_models;       % total number of cells sampled
    case 'GP'
        load([control_path control_list{1,2}{1}],'mean_GPe','GPe_ISIhist')
        ISI_hist_length = length(GPe_ISIhist{1});
        n_control = length(mean_GPe) * n_control_models;       % total number of cells sampled
        
        load([d2_path d2_list{1,2}{1}],'mean_GPe')
        n_d2 = length(mean_GPe) * n_d2_models;       % total number of cells sampled
    otherwise
        error('Unknown analysis requested');
        
end

n_loop = min(n_control_batches,n_d2_batches); 

% storage
n_control_gamma_peaks = zeros(n_loop,1);
n_d2_gamma_peaks = zeros(n_loop,1);
ISI_diff = zeros(n_loop,ISI_hist_length);

% loop through batches and compute...
cl_fig
for loop1 = 1:n_loop
    da03_firing = zeros(n_control,1); 
    da03_peaks = [];
    da03_isi = zeros(ISI_hist_length,1);
    da03_peak_hist = zeros(length(x),1);
    da03_firing_hist = zeros(length(x),1);

    da1_firing = zeros(n_d2,1); 
    da1_peaks = [];
    da1_isi = zeros(ISI_hist_length,1);
    da1_peak_hist = zeros(length(x),1);
    da1_firing_hist = zeros(length(x),1);
    
    % load control combined analysis
    load([control_path control_list{loop1,1}]) 
    
    switch analyse
        case 'STN'
            isi_rates = isi_STN_rates;
            isi_hist =  isi_STN_hist;
            batch_sig_fs_list = STN_sig_fs_list;
        case 'GP'
            isi_rates = isi_GPe_rates;
            isi_hist =  isi_GPe_hist;
            batch_sig_fs_list = GPe_sig_fs_list;
    end
    
    for loop2 = 1:n_control_models
        da03_peaks = [da03_peaks batch_sig_fs_list{loop2}];
    end
    
    for loop2 = 1:n_control
        da03_firing(loop2) = mean(isi_rates{loop2}); 
        da03_isi = da03_isi + isi_hist{loop2}./n_control; 
    end

    n_control_gamma_peaks(loop1) = length(find(da03_peaks >= 40 & da03_peaks <= 80)); 
    da03_peak_hist = histc(da03_peaks,x);   % histogram of all sig peak freqs
    da03_firing_hist = histc(da03_firing,x);

     % load D2 combined analysis
    load([d2_path d2_list{loop1,1}]) 
    
    switch analyse
        case 'STN'
            isi_rates = isi_STN_rates;
            isi_hist =  isi_STN_hist;
            batch_sig_fs_list = STN_sig_fs_list;
        case 'GP'
            isi_rates = isi_GPe_rates;
            isi_hist =  isi_GPe_hist;
            batch_sig_fs_list = GPe_sig_fs_list;
    end

    for loop2 = 1:n_d2_models
        da1_peaks = [da1_peaks batch_sig_fs_list{loop2}];
    end
    
    for loop2 = 1:n_d2
        da1_firing(loop2) = mean(isi_rates{loop2}); 
        da1_isi = da1_isi + isi_hist{loop2}./n_d2; % mean ISI hist
    end
    
    n_d2_gamma_peaks(loop1) = length(find(da1_peaks >= 40 & da1_peaks <= 80)); 
    da1_peak_hist = histc(da1_peaks,x);   % histogram of all sig peak freqs
    da1_firing_hist = histc(da1_firing,x);     % histogram of mean firing rates

    % find change in ISI values...
    ISI_diff(loop1,:) = ((da1_isi ./ da03_isi) .* 100) - 100;    % percentage change in mean ISI    
    
    
    % plot stuff
%     figure(loop1)
%     subplot(711),bar(x,STN_da03_peak_hist,'histc');
%     subplot(712),bar(x,STN_da1_peak_hist,'histc');
% 
%     subplot(713),bar(x,STN_da03_firing_hist,'histc');
%     subplot(714),bar(x,STN_da1_firing_hist,'histc');
% 
%     subplot(715),bar(STN_x_hist,STN_da03_isi);
%     subplot(716),bar(STN_x_hist,STN_da1_isi);
%     
%     subplot(717),bar(STN_x_hist,ISI_diff);
    
    x = x';
    x_hist = STN_x_hist';
    temp_diff = ISI_diff(loop1,:)';
    save(['batch_' num2str(loop1) '_Gamma_x.txt'],'x','-ascii')
    save(['batch_' num2str(loop1) '_Gamma_x_hist.txt'],'x_hist','-ascii')
    save(['batch_' num2str(loop1) '_' analyse '_Gamma_ISI_diff.txt'],'temp_diff','-ascii')
    save(['batch_' num2str(loop1) '_' analyse '_da03_peak_hist.txt'],'da03_peak_hist','-ascii')
    save(['batch_' num2str(loop1) '_' analyse '_da1_peak_hist.txt'],'da1_peak_hist','-ascii')
    save(['batch_' num2str(loop1) '_' analyse '_da03_firing_hist.txt'],'da03_firing_hist','-ascii')
    save(['batch_' num2str(loop1) '_' analyse '_da1_firing_hist.txt'],'da1_firing_hist','-ascii')
    save(['batch_' num2str(loop1) '_' analyse '_da03_isi.txt'],'da03_isi','-ascii')
    save(['batch_' num2str(loop1) '_' analyse '_da1_isi.txt'],'da1_isi','-ascii')


end %% batches loop

mean_ISI_diff = mean(ISI_diff)';
std_ISI_diff = std(ISI_diff)';

save([analyse '_Gamma_n_control_peaks.txt'],'n_control_gamma_peaks','-ascii');
save([analyse '_Gamma_n_d2_peaks.txt'],'n_d2_gamma_peaks','-ascii');
save([analyse '_Gamma_all_ISI_diff.txt'],'ISI_diff','-ascii');
save([analyse '_Gamma_mean_ISI_diff.txt'],'mean_ISI_diff','-ascii');
save([analyse '_Gamma_std_ISI_diff.txt'],'std_ISI_diff','-ascii');
