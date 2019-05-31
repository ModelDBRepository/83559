function c_fname = combine_stngpe_analysis(file_list,flags_fname,pathroot,varargin)

% COMBINE_STNGPE_ANALYSIS summarise results of model batch
%   C = COMBINE_STNGPE_ANALYSIS(L,F,PATH,P) combines the results from every
%   individual analysis file in list L (a cell array), according to the analysis flags set
%   in file F (a string). Results are saved to a file in path PATH,
%   uniquely named and with optional prefix P.

prefix = '';
if nargin >= 4 
    if ~isstr(varargin{1})
        error('Supplied filename prefix is not a string')
    else
        prefix = varargin{1}; 
    end
end

eval(flags_fname)   % evaluate analysis flag file

% generate file name
time_now = clock;
unique_name = datestr(time_now,30);  

c_fname = [prefix '_' unique_name '_combined'];
path_fname = [pathroot c_fname '.mat'];

% variables 
n_models = length(file_list);

% storage
mean_STN_list = [];
mean_GPe_list = [];
mean_GPi_list = [];

isi_STN_rates = {};
isi_STN_hist = {};
isi_STN_x = {};
isi_GPe_rates = {};
isi_GPe_hist = {};
isi_GPe_x = {};

total_STN_bursts = 0;
total_GPe_bursts = 0;

STN_powers = {};
STN_freqs = {};
STN_LFOs = {};
STN_sig_fs_list = {};

GPe_powers = {};
GPe_freqs = {};
GPe_LFOs = {};
GPe_sig_fs_list = {};

GPe_avg_wvfrm_list = [];
STN_avg_wvfrm_list = [];
cortex_mean_list = [];

%%%% combine data %%%%%%%

% exract from files and compile
for loop = 1:n_models
    
    load([pathroot file_list{loop}]);
    
    if do_means
        mean_STN_list = [mean_STN_list; mean_STN'];
        mean_GPe_list = [mean_GPe_list; mean_GPe'];
        mean_GPi_list = [mean_GPi_list; mean_GPi'];       
    end
    
    if do_isi
        isi_STN_rates = [isi_STN_rates; STN_rates];
        isi_STN_hist = [isi_STN_hist; STN_ISIhist];
        isi_STN_x = [isi_STN_x; STN_x_isi];
        
        isi_GPe_rates = [isi_GPe_rates; GPe_rates];
        isi_GPe_hist = [isi_GPe_hist; GPe_ISIhist];
        isi_GPe_x = [isi_GPe_x; GPe_x_isi];
    end

    if do_spectra
        STN_powers = [STN_powers {STN_spect_res.powers}];
        STN_freqs = [STN_freqs {STN_spect_res.freqs}];
        STN_LFOs = [STN_LFOs {STN_spect_res.LFO}];
        STN_sig_fs_list = [STN_sig_fs_list {STN_sig_fs}];
        % create storage for sig freqs...
        
        GPe_powers = [GPe_powers {GPe_spect_res.powers}];
        GPe_freqs = [GPe_freqs {GPe_spect_res.freqs}];
        GPe_LFOs = [GPe_LFOs {GPe_spect_res.LFO}];
        GPe_sig_fs_list = [GPe_sig_fs_list {GPe_sig_fs}];
    end
    
    if do_spike_trig_avg
        % has already been computed per "animal" so cannot be further
        % analysed - just combine results
        STN_avg_wvfrm_list = [STN_avg_wvfrm_list; STN_avg_wvfrm];
        GPe_avg_wvfrm_list = [GPe_avg_wvfrm_list; GPe_avg_wvfrm];
        cortex_mean_list = [cortex_mean_list; ctx_mean];
    end
end

%% combine analyses %%%%
n_STN = length(mean_STN) * n_models;    % total number of STN cells sampled
n_GPe = length(mean_GPe) * n_models;    % total number of GPe cells sampled
n_GPi = length(mean_GPi) * n_models;    % total number of GPi cells sampled

if do_means
    % find mean across neurons in each nucleus
    STN_Hz = mean(mean_STN_list);
    GPe_Hz = mean(mean_GPe_list);
    GPi_Hz = mean(mean_GPi_list);         % included for tonic firing rate setting

    % find std errors of firing rates in each nucleus
    std_STN = std(mean_STN_list);
    std_GPe = std(mean_GPe_list);
    std_GPi = std(mean_GPi_list);

    sem_STN = std_STN / sqrt(n_STN);
    sem_GPe = std_GPe / sqrt(n_GPe);
    sem_GPi = std_GPi / sqrt(n_GPi) ;       % computes the standard error of the mean (SEM) as quoted in e.g Urbain et al 2000

    % find CV
    CV_STN = std_STN ./ STN_Hz;
    CV_GPe = std_GPe ./ GPe_Hz;
    CV_GPi = std_GPi ./ GPi_Hz;

    fprintf(1, 'mean STN rate %.2f:  stderr of mean STN rate %.2f: std dev of STN rate %.2f\n', STN_Hz, sem_STN, std_STN);
    fprintf(1, 'mean GPe rate %.2f:  stderr of mean GPe rate %.2f: std dev of STN rate %.2f\n', GPe_Hz, sem_GPe, std_GPe);
    fprintf(1, 'mean GPi rate %.2f:  stderr of mean GPi rate %.2f: std dev of STN rate %.2f\n', GPi_Hz, sem_GPi, std_GPi);
end

if do_isi
    
end

if do_spectra
    % create mean of all periodograms
    N_freqs = cellfun('length', STN_freqs);
    N_non_null_spectsSTN = length(N_freqs);
    [max_N_fs, index] = max(N_freqs);
    STN_mu_power = zeros(max_N_fs, 1);
    STN_LFO_count = 0;

    STN_freqs_base = STN_freqs{index};
    for loop = 1:N_non_null_spectsSTN
        if N_freqs(loop) > 0
            if N_freqs(loop) ~= max_N_fs
                powers = interp1(STN_freqs{loop}, STN_powers{loop}, STN_freqs_base);
                powers = powers';
            else
                powers = STN_powers{loop};
            end
        if STN_LFOs{loop}
           STN_LFO_count = STN_LFO_count + 1;
        end 
        else
            powers = zeros(max_N_fs, 1);
        end
        STN_mu_power = STN_mu_power + powers;
    end
    STN_mu_power = STN_mu_power ./ N_non_null_spectsSTN;
    
    N_freqs = cellfun('length', GPe_freqs);
    N_non_null_spectsGP = length(N_freqs);
    [max_N_fs, index] = max(N_freqs);
    GPe_mu_power = zeros(max_N_fs,1);
    GPe_LFO_count = 0;

    GPe_freqs_base = GPe_freqs{index};
    for loop=1:N_non_null_spectsGP
        if N_freqs(loop) > 0
            if N_freqs(loop) ~= max_N_fs
                powers = interp1(GPe_freqs{loop}, GPe_powers{loop}, GPe_freqs_base);
                powers = powers';
            else
                powers = GPe_powers{loop};
            end
        if GPe_LFOs{loop}
           GPe_LFO_count = GPe_LFO_count + 1;
        end 
        else
            powers = zeros(max_N_fs, 1);
        end
        GPe_mu_power = GPe_mu_power + powers;
    end
    GPe_mu_power = GPe_mu_power ./ N_non_null_spectsGP;
    
    if do_display
        figure
        subplot(2,1,1)
        plot(STN_freqs_base, STN_mu_power);
        title('Mean STN periodogram')
        subplot(2,1,2)
        plot(GPe_freqs_base, GPe_mu_power);
        title('Mean GP periodogram')
        
        drawnow
    end
    
end
%%% save combined data
save(path_fname,'file_list','n_STN','n_GPe','n_GPi');

if do_means
    save(path_fname, 'STN_Hz','GPe_Hz','GPi_Hz',...
    'sem_STN', 'std_STN', 'sem_GPe', 'std_GPe', 'sem_GPi', 'std_GPi', 'CV_STN', 'CV_GPe', 'CV_GPi', '-append');
end

if do_isi
    save(path_fname, 'isi_STN_rates','isi_STN_hist','isi_STN_x','STN_x_hist',...
        'isi_GPe_rates','isi_GPe_hist','isi_GPe_x','GPe_x_hist', '-append');
end

if do_spectra
    save(path_fname, 'STN_powers', 'GPe_powers','STN_freqs','GPe_freqs',...
        'STN_freqs_base', 'STN_mu_power', 'GPe_freqs_base', 'GPe_mu_power', 'GPe_LFOs', ...
        'GPe_LFO_count', 'STN_LFOs', 'STN_LFO_count','STN_sig_fs_list', 'GPe_sig_fs_list','-append');
end

if do_spike_trig_avg
    save(path_fname, 'cortex_mean_list', 'GPe_avg_wvfrm_list','STN_avg_wvfrm_list', '-append');
end
