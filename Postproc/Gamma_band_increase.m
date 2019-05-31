%% script to compute Gamma band increase using data from models-as-animals
%% protocol
%% Mark Humphries 11/5/2006

clear all

control_batch = 'Gamma_a_20060407T204507_batch.mat';
%d2_batch = 'Gamma_b_20060407T152129_batch.mat';
% treat NMDA manipulation as equivalent to D2 agonist
d2_batch = 'Gamma_noNMDAinGP_a_20060512T211003_batch.mat';

control_path = '../ResultsArchive/Gamma-band/CompleteModel/ConditionA/';
% d2_path = '../ResultsArchive/Gamma-band/CompleteModel/ConditionB/';
% treat NMDA manipulation as equivalent to D2 agonist
d2_path = '../ResultsArchive/Gamma-band/NoNMDAinGP/';

analyse = 'STN'; % either 'STN' or 'GP';


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
load([control_path control_batch]);
control_list = batch_analysis_list;
[n_control_batches c] = size(control_list);
n_control_models = length(control_list{1,2});

load([d2_path d2_batch]);
d2_list = batch_analysis_list;
[n_d2_batches c] = size(d2_list);
n_d2_models = length(d2_list{1,2});

% set number of loops to traverse - number of batches
n_loop = min(n_control_batches,n_d2_batches); 

n_model_loop = min(n_control_models,n_d2_models);   % n_d2_models should be less if followed Brown et al

gamma_mean_ratio = zeros(n_loop,1);
gamma_se_ratio = zeros(n_loop,1);


% load first combined analysis to get size of spectra....
load([control_path control_list{1,1}])
n_freqs = length(STN_mu_power);
mean_batch_control_spectra = zeros(n_loop,n_freqs);
mean_batch_d2_spectra = zeros(n_loop,n_freqs);

cl_fig
for loop1 = 1:n_loop

    % load control combined analysis
    load([control_path control_list{loop1,1}])
    
    switch analyse
        case 'STN'
            n_cells_control_model = n_STN / n_control_models;
            freqs_base1 = STN_freqs_base;
            mean_batch_control_spectra(loop1,:) = STN_mu_power;
            powers = STN_powers; 
            freqs = STN_freqs;
        case 'GP'
            n_cells_control_model = n_GPe / n_control_models;
            freqs_base1 = GPe_freqs_base;
            mean_batch_control_spectra(loop1,:) = GPe_mu_power;
            powers = GPe_powers; 
            freqs = GPe_freqs;
        otherwise
            error('unknown analysis selection');
    end
    
    % loop through powers, computing means-per-model
    mean_control_model_spectra = zeros(n_model_loop,n_freqs);
    for loop2 = 1:n_model_loop
        for loop3 = 1:n_cells_control_model
            idx = (loop2-1) * n_cells_control_model + loop3;
            if powers{idx}
                if length(powers{idx}) ~= n_freqs   % if for some reason there are less frequencies than expected
                    this_powers = interp1(freqs{idx}, powers{idx}, freqs_base1);
                    % powers = powers;
                else
                    this_powers = powers{idx}';
                end
                mean_control_model_spectra(loop2,:) = mean_control_model_spectra(loop2,:) + this_powers ./ n_cells_control_model;
            end
        end
    end

    % load D2 combined analysis
    load([d2_path d2_list{loop1,1}]) 
    
    switch analyse
        case 'STN'
            n_cells_d2_model = n_STN / n_d2_models;
            freqs_base2 = STN_freqs_base;
            mean_batch_d2_spectra(loop1,:) = STN_mu_power;
            powers = STN_powers; 
            freqs = STN_freqs;
        case 'GP'
            n_cells_d2_model = n_GPe / n_d2_models;
            freqs_base2 = GPe_freqs_base;
            mean_batch_d2_spectra(loop1,:) = GPe_mu_power;
            powers = GPe_powers; 
            freqs = GPe_freqs;
        otherwise
            error('unknown analysis selection');
    end
    
    % loop through powers, computing means-per-model
    mean_d2_model_spectra = zeros(n_model_loop,n_freqs);
    for loop2 = 1:n_model_loop
        for loop3 = 1:n_cells_d2_model
            idx = (loop2-1) * n_cells_d2_model + loop3;
            if powers{idx}
                if length(powers{idx}) ~= n_freqs   % if for some reason there are less frequencies than expected
                    this_powers = interp1(freqs{idx}, powers{idx}, freqs_base2);
                    %powers = powers;
                else
                    this_powers = powers{idx}';
                end

                mean_d2_model_spectra(loop2,:) = mean_d2_model_spectra(loop2,:) + this_powers ./ n_cells_d2_model;
            end
        end
    end

    %max_power = max([max(STN_power1) max(STN_power2)]);
    
    
%     figure(loop1)
%     h1 = plot(STN_freqs1,STN_power1,'k:','LineWidth',1);
%     hold on
%     h2 = plot(STN_freqs2,STN_power2,'k','LineWidth',1);
%     xlabel('Frequency (Hz)')
%     ylabel('Power')
%     axis([0 100 0 max_power+5])
% 
%     % alter fonts etc
%     ph1 = get(h1,'Parent');
%     fh1 = get(ph1,'Parent');
%     xlbl1 = get(ph1,'XLabel');
%     ylbl1 = get(ph1,'YLabel');
% 
%     set(ph1,'FontSize',14)
%     set(xlbl1,'FontSize',14);
%     set(ylbl1,'FontSize',14);

    % ph2 = get(h2,'Parent');
    % fh2 = get(ph2,'Parent');
    % set(ph2,'FontSize',14)

    % get all indexes in Gamma range
    gamma = find(freqs_base1 >= 40 & freqs_base1 <= 80);

    % total power in each model's spectra in that range
    gamma_pwr1 = sum(mean_control_model_spectra(:,gamma)');
    gamma_pwr2 = sum(mean_d2_model_spectra(:,gamma)');

    % ratio
    gamma_ratio = gamma_pwr2 ./ gamma_pwr1;
    gamma_mean_ratio(loop1) = mean(gamma_ratio);
    gamma_se_ratio(loop1) = std(gamma_ratio) / sqrt(n_model_loop);
    
    % print graphs
    %print(fh1,'-dpng','-r600','gamma_da02.png');
    %print(fh2,'-dpng','-r600','gamma_da08.png');
    freqs_base1 = freqs_base1';
    freqs_base2 = freqs_base2';
    
    
    mean_d2_model_spectra = mean_d2_model_spectra';
    mean_control_model_spectra = mean_control_model_spectra';
    
    save(['batch_' num2str(loop1) '_' analyse '_freqs_control.txt'],'freqs_base1','-ascii');
    save(['batch_' num2str(loop1) '_' analyse '_freqs_d2.txt'],'freqs_base2','-ascii');
    save(['batch_' num2str(loop1) '_' analyse '_mean_control_model_spectra.txt'],'mean_control_model_spectra','-ascii');
    save(['batch_' num2str(loop1) '_' analyse '_mean_d2_model_spectra.txt'],'mean_d2_model_spectra','-ascii');

end

overall_mean_control_spectra = mean(mean_batch_control_spectra);
overall_mean_control_spectra_std = std(mean_batch_control_spectra);
overall_mean_d2_spectra = mean(mean_batch_d2_spectra);
overall_mean_d2_spectra_std = std(mean_batch_d2_spectra);

overall_mean_control_spectra = overall_mean_control_spectra';
overall_mean_control_spectra_std = overall_mean_control_spectra_std';
overall_mean_d2_spectra = overall_mean_d2_spectra';
overall_mean_d2_spectra_std = overall_mean_d2_spectra_std';

mean_batch_control_spectra = mean_batch_control_spectra';
mean_batch_d2_spectra = mean_batch_d2_spectra';

save([analyse '_overall_mean_control_spectra.txt'],'overall_mean_control_spectra','-ascii');
save([analyse '_overall_mean_control_spectra_std.txt'],'overall_mean_control_spectra_std','-ascii');
save([analyse '_overall_mean_d2_spectra.txt'],'overall_mean_d2_spectra','-ascii');
save([analyse '_overall_mean_d2_spectra_std.txt'],'overall_mean_d2_spectra_std','-ascii');
%save([analyse '_mean_batch_control_spectra.txt'], 'mean_batch_control_spectra','-ascii');
%save([analyse '_mean_batch_d2_spectra.txt'], 'mean_batch_d2_spectra','-ascii');
save([analyse '_gamma_mean_ratio.txt'], 'gamma_mean_ratio','-ascii');
save([analyse '_gamma_se_ratio.txt'], 'gamma_se_ratio','-ascii');



%%% do stuff here to look at spread over repeated experiments (batches)




