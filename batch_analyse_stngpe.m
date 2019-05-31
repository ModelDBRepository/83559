function batch_analyse_stngpe(e_fname,flags_fname,pathroot,an_fname)

% BATCH_ANALYSE_STNGPE to analyse firing statistics of sampled STN and GPe neurons
%   BATCH_ANALYSE_STNGPE(E,F,A) analyses the extracted data in file E (a
%   string) using the analysis methods specifed in flags file F (a string).
%   Saves the results to the file named A.
%
% Computes in order for STN and GPe:
% (0) Plots raster and single-neuron trace, if requested
% (1) Mean, SD, SE, and CV of firing rates over whole simulation (includes SNr)
% (2) ISIs, ISI histograms (including input)
% (3) Burst firing analysis (KBSTA method)
% (4) Power-spectra (multi-taper method - chronux toolbox)
%     (a) finds most significant peak
%     (b) finds any significant peak in LFO (<1.5Hz) range (Magill et al. 2001)
% (5) Intra- and inter-nucleus cross-correlation, using peak/trough detection method set at start
% (6) Auto-correlograms
% (7) Coherency (chronux toolbox)
% (8) Spike-triggered averaging (including computation of pseudo-EEG on which it's based)
%
% Mark Humphries & Kevin Gurney 6/1/2006. 

load([pathroot e_fname])       % load extracted results

eval(flags_fname)   % evaluate analysis flag file

 [a s] =  system('hostname');
 host = findstr(s, 'iceberg');
 if ~isempty(host) % on iceberg
     ice_path1 = genpath('/home1/pc/pc1mdh/BG spiking model');
     ice_path2 = genpath('/home1/pc/pc1mdh/Matlab Tools');
     path(path, ice_path1);
     path(path, ice_path2);
 end

cl_fig

do_display = 0; % set for showing graphs interactively

% t_out_t = double(out_t) .* dt; % timestamps for network events
t_in_t = double(in_t) .* dt;    % timestamps for input streams
time_steps = time_seconds / dt;

clear link_n link_wAMPA link_wNMDA link_wGABAa

%%%%%%%%%%%%%%%% SET ANALYSIS PARAMETERS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for ISI data
ISIbins = 100; %No. of bins in ISI hist.
limits = [0.005 0.6];  % ISI hist limits in seconds (outside these bounds get put into end bins)

% parameters for burst analysis
eta_on = 0.4;
eta_off = 0.4;

% parameters for periodogram
freq_range = [0.1 100];    % min and max freqs in periodogram
tapers = [3 5];                         % recommended in Chronuz; 5 tapers with bandwidth*time = 3
pad = 2;                                % Padding the frequency sampling by another 2 powers of 2;
Fs = 1 ./ dt;                           % best to use the underlying sampling;
err_bars = [1 0.01];                    % method 1 (theoretical) and p_sig = 0.01
fscorr = 0;                             %don't use finite size corrections
Min_spikes = 12;                        % have to have at least this number of spikes to do spectrum

LFO_freq = 1.5; % (For Magill's analysis)

% parameters for auto- and cross-correlogram
binsize = 0.001;   % in seconds
maxlag = 0.5;       % in seconds
acf_maxlag = 2;     % window size for auto-correlograms

p = 0.01;           % 99% confidence interval
num_bins = 2*maxlag/binsize + 1;
method1 = 0; % method 1 (proportion of bins outside confidence limits) or method 2 (Raz et al 2000)
if method1
    num_sig_bins = round(num_bins * p); % number of bins over confidence interval required - 
    % use for testing any bins
else
    num_sig_bins = 3;           % when using consecutive bins, more than 3 is significant?                                    
end

% parameters for coherence
freq_rangeC = [0.1 100];    % min and max freqs in periodogram
tapersC = [3 5];                         % recommended in Chronuz; 5 tapers with bandwidth*time = 3
padC = 2;                                % Padding the frequency sampling by another 2 powers of 2;
FsC = 1 ./ dt;                           % best to use the underlying sampling;
err_barsC = [2 0.01];                    % method (jacknife) and p_sig = 0.01
fscorrC = 0;                             %don't use finite size corrections
Min_spikesC = 12;                        % have to have at least this number of spikes to do spectrum


%%% parameters for spike-triggered averaging
    
% a. parameters for smoothed input firing rate estimation
    smooth_win_size = 0.01;     % width in seconds
    step_size = 10;             % window step in time-steps
    window = 'alpha';           % alt: 'gaussian'
% b. parameters for spike-triggered average
    trig_win_half = 1;          % half-size of averaging window in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% some basic pre-processing %%%%%%
input_times = cell(neurons_per_nucleus,1);

% create time-stamp arrays for each input
for loop=1:neurons_per_nucleus    
    input_times{loop} = dt .* double(in_t(in_n == EXT(loop)));
end

%%%% create data structures for periodogram results - remove this later if
%%%% we tidy up data storage!

% number of neurons from each structure
n_STN = length(STN_times);
n_GPe = length(GPe_times);
n_GPi = nan;

for loop = 1:n_STN
    STN_data(loop) = struct('times', STN_times{loop}');     % for multi_taper - if multichannel use is made of the
                                                            % chronux call, then need this format for data
end

for loop = 1:n_GPe
    GPe_data(loop) = struct('times', GPe_times{loop}');   % for multi_taper
end

%%%%%%%%%% RASTER PLOT %%%%%%%%%%%%%%%%%%
if do_display & do_raster
    structured_raster(in_n, in_t, out_n, out_t, dt, neurons_per_nucleus, SD1, SD2, STN, GPe, GPi, time_seconds);
    
    display('press any key to continue');
    pause;
end

%%%%%%%%%%%% trace %%%%%%%%%%%%%%%%%%%
% show trace if appropriate
if any([trace_n == STN trace_n == GPe]) & do_trace & do_display
	trace_times = 1:1:time_steps;
	figure
	plot(trace_times,trace_vals(:,1))
	figure
	plot(trace_times,trace_vals(:,2))
	figure
	plot(trace_times,trace_vals(:,3))
    tile
    drawnow
end

%%%%%%%% basic mean firing rate stats
if do_means
    
    % determine mean outputs (firing rates) for each neuron
    mean_STN = cellfun('length',STN_times) ./ time_seconds;     % length of each array is number of spikes! 
    mean_GPe = cellfun('length',GPe_times) ./ time_seconds;
   
    
    % find mean across neurons in each nucleus
    STN_Hz = mean(mean_STN);
    GPe_Hz = mean(mean_GPe);
    
    % find std errors of firing rates in each nucleus
    std_STN = std(mean_STN);
    std_GPe = std(mean_GPe);
    
    
    sem_STN = std_STN / sqrt(neurons_per_nucleus);
    sem_GPe = std_GPe / sqrt(neurons_per_nucleus);
    
    % find CV of firing rates for each nucleus (*not* neuron)
    % for analysis of single neuron trains, find CV of ISIs...
    CV_STN = std_STN ./ STN_Hz;
    CV_GPe = std_GPe ./ GPe_Hz;
    
    if exist('GPi_times')
        n_GPi = length(GPi_times);
        mean_GPi = cellfun('length',GPi_times) ./ time_seconds;
        GPi_Hz = mean(mean_GPi);         % included for tonic firing rate setting
        std_GPi = std(mean_GPi);
        sem_GPi = std_GPi / sqrt(neurons_per_nucleus) ;       % computes the standard error of the mean (SEM) as quoted in e.g Urbain et al 2000
        CV_GPi = std_GPi ./ GPi_Hz;
    else
        n_GPi = nan;
        mean_GPi = nan;
        GPi_Hz = nan;
        std_GPi = nan;
        sem_GPi = nan;
        CV_GPi = nan;
    end 
    
%     fprintf(1, 'mean STN rate %.2f:  stderr of mean STN rate %.2f: std dev of STN rate %.2f\n', STN_Hz, sem_STN, std_STN);
%     fprintf(1, 'mean GPe rate %.2f:  stderr of mean GPe rate %.2f: std dev of STN rate %.2f\n', GPe_Hz, sem_GPe, std_GPe);
%     fprintf(1, 'mean GPi rate %.2f:  stderr of mean GPi rate %.2f: std dev of STN rate %.2f\n', GPi_Hz, sem_GPi, std_GPi);

end
if do_display
    display('press any key to continue 2');
    pause;
end
 


%%%%%%%%%%%%%%%% isi %%%%%%%%%%%%%%%%%%%%%
STN_rates = cell(n_STN,1);            % instantaneous rates for *all* events (not just bursts)
STN_x_isi = cell(n_STN,1);            % time stamps for rates above (almost just timestamps )
STN_ISIhist = cell(n_STN,1);          % ISI histograms

GPe_rates = cell(n_GPe,1);
GPe_x_isi = cell(n_GPe,1);
GPe_ISIhist = cell(n_GPe,1);

input_rates = cell(length(EXT),1);
input_x_isi = cell(length(EXT),1);
input_ISIhist = cell(length(EXT),1);

cum_STN_spikes =  sparse(1,time_steps);
cum_GPe_spikes =  sparse(1,time_steps);
cum_input_spikes =  sparse(1,time_steps);
if do_isi
    for loop = 1:n_STN
        % generate ISI values
        [STN_rates{loop},STN_ISIhist{loop},STN_x_isi{loop},STN_x_hist] = LIF_ISI_analysis(STN_times{loop},ISIbins,limits);
        STN_discrete_times = sparse(1,time_steps);
        STN_discrete_times(round(STN_times{loop} ./ dt)) = 1;
        cum_STN_spikes = cum_STN_spikes + STN_discrete_times;
    end

    for loop = 1:n_GPe
        [GPe_rates{loop},GPe_ISIhist{loop},GPe_x_isi{loop},GPe_x_hist] = LIF_ISI_analysis(GPe_times{loop},ISIbins,limits);
        GPe_discrete_times = sparse(1,time_steps);
        GPe_discrete_times(round(GPe_times{loop} ./ dt)) = 1;
        cum_GPe_spikes = cum_GPe_spikes + GPe_discrete_times;
    end
    
    for loop = 1:neurons_per_nucleus
        [input_rates{loop},input_ISIhist{loop},input_x_isi{loop},input_x_hist] = LIF_ISI_analysis(input_times{loop},ISIbins,limits);
        input_discrete_times = sparse(1,time_steps);
        input_discrete_times(round(input_times{loop} ./ dt) + 1) = 1; % input time can be zero causes a problem
        input_discrete_times = input_discrete_times(1:time_steps);
        cum_input_spikes = cum_input_spikes + input_discrete_times;
    end
    rs_STN = find_mean_isis(cum_STN_spikes);
    rs_GPe = find_mean_isis(cum_GPe_spikes);
    rs_in = find_mean_isis(cum_input_spikes);
    
    
    %%%% show individual isi pots for STN %%%%%
    No_subplots_isi = 4;
    if do_display
        done = 0;
        loop = 0;
        while loop < n_STN
            plot_number = mod(loop, No_subplots_isi) + 1;
            if  plot_number == 1
                figure
            end
            subplot(No_subplots_isi, 1, plot_number);
            plot(STN_x_isi{loop+1}, STN_rates{loop+1});
            drawnow;
            loop = loop + 1;
        end
    end
    
    if do_display
       display('press any key to continue');
       pause
    end
    %%%%%%%%%%%%%%%%%%%%
    
    
    
    if do_display
        figure
        plot(rs_STN);
        title('mean isis for STN')
        figure
        plot(rs_in);
        title('mean isis for input')
    end
end



%%%%%%%%%%%%%% bursts %%%%%%%%%%%%%%%
if do_bursts
    % determine burst status
    STN_bursts = cell(n_STN,1);           % burst start and end times
    STN_burst_isis = cell(n_STN,1);       % within burst ISIs (entire sequence for each burst)
    
    
    % repeat for GPe....
    GPe_bursts = cell(n_GPe,1);           % burst start and end times
    GPe_burst_isis = cell(n_GPe,1);       % within burst ISIs
    
    
    
    for loop = 1:n_STN
        % burst-detection
        fprintf(1, 'STN burst loop count %d\n', loop);
        
        [STN_bursts{loop},STN_burst_isis{loop}] = kbsta(STN_times{loop},[],eta_on,eta_off);
        [N_STN_bursts(loop) c] = size(STN_bursts{loop});
    end
    
    for loop = 1:n_GPe
        % burst-detection
        fprintf(1, 'GPe burst loop count %d\n', loop);
       
        [GPe_bursts{loop},GPe_burst_isis{loop}] = kbsta(GPe_times{loop},[],eta_on,eta_off);
        [N_GPe_bursts(loop) c] = size(GPe_bursts{loop});        
    end
    
    if do_display
        figure
        subplot(2,1,1)
        hist(N_STN_bursts, [0:15]);
        title('histogram of STN bursts');
        ylabel('No of cells');
        xlabel('No bursts');
        
        subplot(2,1,2)
        hist(N_GPe_bursts, [0:15]);
        title('histogram of GPe bursts');
        ylabel('No of cells');
        xlabel('No bursts');
        
        drawnow
    end

end

%%%%%%%%%%% spectra %%%%%%%%%%%%%%%%%%

% periodogram analysis
% [STN_results,STN_store_power,STN_fs_array] = scargle_analysis(STN_rates,STN_x_isi,dt, time_seconds, [], freqs);
% STNnans = isnan(STN_results(:,1)) | STN_results(:,1)==0;
% STN_sum_results = STN_results;
% STN_sum_results(STNnans,:) = [];       % results matrix without empty records
% 
% [GPe_results,GPe_store_power,GPe_fs_array] = scargle_analysis(GPe_rates,GPe_x_isi,dt, time_seconds, [],freqs);
% GPenans = isnan(GPe_results(:,1)) | GPe_results(:,1)==0;
% GPe_sum_results = GPe_results;
% GPe_sum_results(GPenans,:) = [];       % results matrix without empty records

if do_display
    display('press any key to continue 3');
    pause;
end
if do_spectra
    % multi-taper analysis
    
    STN_spect_res = struct('freqs', [], 'powers', [], 'errs', [], 'peak_freq', [], 'peak_power', [], 'LFO', []);
    GPe_spect_res = struct('freqs', [], 'powers', [], 'errs', [], 'peak_freq', [], 'peak_power', [], 'LFO', []);
    
    disp('doing multi-taper spectral analysis');
    
    for loop = 1:n_STN
        fprintf(1, 'STN spectral loop count %d\n', loop);
        
        %% STN%%% 
        data = STN_data(loop).times;
        [r c] = size(data);
        if r >= Min_spikes
            [STN_MT_power, STN_MT_fs_array, R, STN_MT_errs] = mtspectrumpt(data, tapers, pad, Fs, freq_range, err_bars, 0, fscorr); % 0 is no trial average
            STN_err_lo = STN_MT_errs(1,:);
            STN_spect_res(loop).freqs = STN_MT_fs_array;
            STN_spect_res(loop).powers = STN_MT_power;
            STN_spect_res(loop).errs = STN_MT_errs;
            
            sig_test_STN_power = STN_err_lo - mean(STN_MT_power);
            [STN_peak STN_peak_index] = max(sig_test_STN_power);
            
            if STN_peak >0
                STN_spect_res(loop).peak_freq = STN_MT_fs_array(STN_peak_index);
                STN_spect_res(loop).peak_power = STN_peak;
            else
                STN_spect_res(loop).peak_freq = NaN;
                STN_spect_res(loop).peak_power = NaN;
            end

            %%%% extra for LFOs in magill work %%%%%%%%%

            LFO_freq_indices = find(STN_MT_fs_array < LFO_freq);
            STN_LFO_sig_test = sig_test_STN_power(LFO_freq_indices);
            [STN_peak_LFO STN_peak_index_LFO] = max(STN_LFO_sig_test);
            if STN_peak_LFO > 0
                STN_spect_res(loop).LFO = 1;
            else
                STN_spect_res(loop).LFO = 0;
            end

	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            STN_spect_res(loop).freqs = [];
            STN_spect_res(loop).powers = [];
            STN_spect_res(loop).errs = [];
            STN_spect_res(loop).peak_freq = NaN;
            STN_spect_res(loop).peak_power = NaN;
            STN_spect_res(loop).LFO = NaN;
        end
    end
    
    for loop = 1:n_GPe
        fprintf(1, 'GPe spectral loop count %d\n', loop);
        %% GPe %%% 
        data = GPe_data(loop).times;
        [r c] = size(data);
        if r >= Min_spikes
            [GPe_MT_power, GPe_MT_fs_array, R, GPe_MT_errs] = mtspectrumpt(data, tapers, pad, Fs, freq_range, err_bars, 0, fscorr); % 0 is no trial average
            GPe_err_lo = GPe_MT_errs(1,:);
            GPe_spect_res(loop).freqs = GPe_MT_fs_array;
            GPe_spect_res(loop).powers = GPe_MT_power;
            GPe_spect_res(loop).errs = GPe_MT_errs;
            
            sig_test_GPe_power = GPe_err_lo - mean(GPe_MT_power);
            [GPe_peak GPe_peak_index] = max(sig_test_GPe_power);
            
            if GPe_peak >0
                GPe_spect_res(loop).peak_freq = GPe_MT_fs_array(GPe_peak_index);
                GPe_spect_res(loop).peak_power = GPe_peak;
            else
                GPe_spect_res(loop).peak_freq = NaN;
                GPe_spect_res(loop).peak_power = NaN;
            end

            %%%% extra for LFOs in magill work %%%%%%%%%

            LFO_freq_indices = find(GPe_MT_fs_array < LFO_freq);
            GPe_LFO_sig_test = sig_test_GPe_power(LFO_freq_indices);
            [GPe_peak_LFO Gpe_peak_index_LFO] = max(GPe_LFO_sig_test);
            if GPe_peak_LFO > 0
                GPe_spect_res(loop).LFO = 1;
            else
                GPe_spect_res(loop).LFO = 0;
            end

	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        else
            GPe_spect_res(loop).freqs = [];
            GPe_spect_res(loop).powers = [];
            GPe_spect_res(loop).errs = [];
            GPe_spect_res(loop).peak_freq = NaN;
            GPe_spect_res(loop).peak_power = NaN;
            GPe_spect_res(loop).LFO = NaN;
        end
    end
    
    sig_STN_peak_freqs = [STN_spect_res.peak_freq]; 
    STN_sig_fs = sig_STN_peak_freqs(find(sig_STN_peak_freqs > 0));
    lows = STN_sig_fs(STN_sig_fs < 5);
    sig_GPe_peak_freqs = [GPe_spect_res.peak_freq];
    GPe_sig_fs = sig_GPe_peak_freqs(find(sig_GPe_peak_freqs > 0));
    lowg = GPe_sig_fs(GPe_sig_fs < 5);

    if do_display
        figure
       
        subplot(4,1,1)
        hist(STN_sig_fs, 20);
        title('significant maximum  frequencies - STN');
        xlabel('peak freq');
        subplot(4,1,2)
        hist(lows, 15);
        title('significant maximum  frequencies in range 0-5Hz - STN');
        xlabel('peak freq');
        
        subplot(4,1,3)
        hist(GPe_sig_fs, 20);
        title('significant maximum  frequencies - Gpe');
        xlabel('peak freq');
        subplot(4,1,4)
        hist(lowg, 15);
        title('significant maximum  frequencies in range 0-5Hz - GPe');
        xlabel('peak freq');
        
        drawnow
    end    
end % do_spectra

%%%%%%%%%%%%%% auto- and cross-correlograms  %%%%%%%%%%
n_STN_pairs = n_STN ^ 2;
n_GPe_pairs = n_GPe ^ 2;
n_inter_pairs = n_STN * n_GPe;

STN_pairs = list_of_pairs(1:n_STN,1:n_STN,'su'); 
GPe_pairs = list_of_pairs(1:n_GPe,1:n_GPe,'su');
STNGPe_pairs = list_of_pairs(1:n_STN,1:n_GPe,'su');

STN_ccf = cell(n_STN_pairs,1);
GPe_ccf = cell(n_GPe_pairs,1);
STNGPe_ccf = cell(n_inter_pairs,1);
total_STN_ccf = zeros(1,num_bins);
total_GPe_ccf = zeros(1,num_bins);
total_STNGPe_ccf = zeros(1,num_bins);
total_STN_pairs = 0;
total_GPe_pairs = 0;
total_STNGPe_pairs = 0;

pos_sigSTNpairs = zeros(n_STN_pairs,1);
neg_sigSTNpairs = zeros(n_STN_pairs,1);
pos_sigGPepairs = zeros(n_GPe_pairs,1);
neg_sigGPepairs = zeros(n_GPe_pairs,1);


if do_xcorr_intra
disp('doing xcorr - intra nucleus');
    % OPTION 2: calculate cross-correlation histograms direct from spike trains

    GPe_peak_and_trough = 0; STN_peak_and_trough = 0; 
    for i = 1:n_STN_pairs
        fprintf(1, 'STN xcorr loop %d\n', i);

        [STN_ccf{i},xSTN,f1,f2,n_pairs] = LIF_xcorr(STN_times{STN_pairs(i,1)},STN_times{STN_pairs(i,2)},binsize,[0 time_seconds],maxlag); 
        total_STN_pairs = total_STN_pairs + n_pairs;
        [mSTN,sSTN,spSTN,stSTN] = xcorr_stats(STN_ccf{i},n_pairs,p);
 
        % test for significance in cross-correlograms
        S_pos_sig_bins = find(STN_ccf{i} >= spSTN);
        S_neg_sig_bins = find(STN_ccf{i} < stSTN);
        
        if method1
            % METHOD 1: number of bins greater than confidence interval exceeds the expected amount
            if length(S_pos_sig_bins) > num_sig_bins     % for just comparing all bins over correlogram
                pos_sigSTNpairs(i) = 1;    
            end
            if length(S_neg_sig_bins) > num_sig_bins 
                neg_sigSTNpairs(i) = 1;   
            end            
        else
            % METHOD 2: number of consecutive significant bins is greater than fixed threshold (3 - Raz et al 2000)
            for loop2 = 1:num_sig_bins
                S_pos_sig_bins = find(diff(S_pos_sig_bins) == 1); % indices of adjacent bins - reduce                
                S_neg_sig_bins = find(diff(S_neg_sig_bins) == 1); % indices of adjacent bins - reduce      
            end
            
            if ~isempty(S_pos_sig_bins)  pos_sigSTNpairs(i) = 1; end   % then there is at least one occurence of a num_sig_bins number of consecutive significant bins
            if ~isempty(S_neg_sig_bins)  neg_sigSTNpairs(i) = 1; end
        end
       
        % omit correlograms with both sig peak and trough 
        if pos_sigSTNpairs(i) == 1 & neg_sigSTNpairs(i) == 1
            pos_sigSTNpairs(i) = 0;
            neg_sigSTNpairs(i) = 0;
            STN_peak_and_trough = STN_peak_and_trough + 1;
        end
        % total CCFs
        total_STN_ccf = total_STN_ccf + STN_ccf{i};

    end
    
    for i = 1:n_GPe_pairs
        fprintf(1, 'GPe xcorr loop %d\n', i);
        
        [GPe_ccf{i},xGPe,f1,f2,n_pairs] = LIF_xcorr(GPe_times{GPe_pairs(i,1)},GPe_times{GPe_pairs(i,2)},binsize,[0 time_seconds],maxlag);
        total_GPe_pairs = total_GPe_pairs + n_pairs;
        [mGPe,sGPe,spGPe,stGPe] = xcorr_stats(GPe_ccf{i},n_pairs,p);
        
        % test for significance in cross-correlograms
        G_pos_sig_bins = find(GPe_ccf{i} >= spGPe);
        G_neg_sig_bins = find(GPe_ccf{i} < stGPe);
        
        if method1
            % METHOD 1: number of bins greater than confidence interval exceeds the expected amount
            if length(G_pos_sig_bins) > num_sig_bins 
                pos_sigGPepairs(i) = 1;   
            end
            if length(G_neg_sig_bins) > num_sig_bins 
                neg_sigGPepairs(i) = 1;   
            end
            
        else
            % METHOD 2: number of consecutive significant bins is greater than fixed threshold (3 - Raz et al 2000)
            for loop2 = 1:num_sig_bins
                G_pos_sig_bins = find(diff(G_pos_sig_bins) == 1); % indices of adjacent bins - reduce                
                G_neg_sig_bins = find(diff(G_neg_sig_bins) == 1); % indices of adjacent bins - reduce                
            end
            
            if ~isempty(G_pos_sig_bins)  pos_sigGPepairs(i) = 1; end 
            if ~isempty(G_neg_sig_bins)  neg_sigGPepairs(i) = 1; end
        end
       
        % omit correlograms with both sig peak and trough 
        if pos_sigGPepairs(i) == 1 & neg_sigGPepairs(i) == 1
            pos_sigGPepairs(i) = 0;
            neg_sigGPepairs(i) = 0;
            GPe_peak_and_trough = GPe_peak_and_trough + 1;
        end

        % total CCFs
        total_GPe_ccf = total_GPe_ccf + GPe_ccf{i};
    end
    
    
    mean_STN_ccf = total_STN_ccf ./ n_STN_pairs;
    mean_GPe_ccf = total_GPe_ccf ./ n_GPe_pairs;
    
    prop_POS_sig_STN_pairs = sum(pos_sigSTNpairs) / n_STN_pairs
    prop_NEG_sig_STN_pairs = sum(neg_sigSTNpairs) / n_STN_pairs
    prop_POS_sig_GPe_pairs = sum(pos_sigGPepairs) / n_GPe_pairs
    prop_NEG_sig_GPe_pairs = sum(neg_sigGPepairs) / n_GPe_pairs
    

    if do_total_xcorr_stats
        [mtotalSTN,stotalSTN,sptotalSTN,sttotalSTN] = xcorr_stats(total_STN_ccf,total_STN_pairs);
        [mtotalGPe,stotalGPe,sptotalGPe,sttotalGPe] = xcorr_stats(total_GPe_ccf,total_GPe_pairs);
        if do_display
            figure
            subplot(211),bar(xSTN,total_STN_ccf,1)
            axis([-0.1 0.1 min(total_STN_ccf)  max(total_STN_ccf)])
            line([-0.1 0.1],[sptotalSTN sptotalSTN],'Color',[1 0 0]);
            subplot(212),bar(xGPe,total_GPe_ccf,1)
            axis([-0.1 0.1 min(total_GPe_ccf)  max(total_GPe_ccf)])
            line([-0.1 0.1],[sptotalGPe sptotalGPe],'Color',[1 0 0]);
            drawnow
        end
    end
end

%%%%%%%%%%%%% do xcorr between nuclei %%%%%%%%%%%%%

pos_sigSTNGPepairs = zeros(n_inter_pairs,1);
neg_sigSTNGPepairs = zeros(n_inter_pairs,1);

if do_xcorr_inter
disp('doing xcorr - inter nucleus');
    STNGPe_peak_and_trough = 0; 

    for i = 1:n_inter_pairs
        [STNGPe_ccf{i},xSTNGPe,f1,f2,n_pairs] = LIF_xcorr(STN_times{STNGPe_pairs(i,1)},GPe_times{STNGPe_pairs(i,2)}, binsize, [0 time_seconds], maxlag);
        total_STNGPe_pairs = total_STNGPe_pairs + n_pairs;
        total_STNGPe_ccf = total_STNGPe_ccf + STNGPe_ccf{i};    
        [mSG,sSG,spSG,stSG] = xcorr_stats(STNGPe_ccf{i},n_pairs,p);

        % test for significance in cross-correlograms
        SG_pos_sig_bins = find(STNGPe_ccf{i} >= spSG);
        SG_neg_sig_bins = find(STNGPe_ccf{i} < stSG);
         if method1
            % METHOD 1: number of bins greater than confidence interval exceeds the expected amount
            if length(SG_pos_sig_bins) > num_sig_bins     % for just comparing all bins over correlogram
                pos_sigSTNGPepairs(i) = 1;    
            end
            if length(SG_neg_sig_bins) > num_sig_bins 
                neg_sigSTNGPepairs(i) = 1;   
            end            
        else
         % METHOD 2: number of consecutive significant bins is greater than fixed threshold (3 - Raz et al 2000)
            for loop2 = 1:num_sig_bins
                SG_pos_sig_bins = find(diff(SG_pos_sig_bins) == 1); % indices of adjacent bins - reduce                
                SG_neg_sig_bins = find(diff(SG_neg_sig_bins) == 1); % indices of adjacent bins - reduce      
                
            end
            
            if ~isempty(SG_pos_sig_bins)  pos_sigSTNGPepairs(i) = 1; end   % then there is at least one occurence of a num_sig_bins number of consecutive significant bins
            if ~isempty(SG_neg_sig_bins)  neg_sigSTNGPepairs(i) = 1; end
        end

        if pos_sigSTNGPepairs(i) == 1 & neg_sigSTNGPepairs(i) == 1
            pos_sigSTNGPepairs(i) = 0;
            neg_sigSTNGPepairs(i) = 0;
            STNGPe_peak_and_trough = STNGPe_peak_and_trough + 1;
        end

    end
    mean_STNGPe_ccf = total_STNGPe_ccf ./ n_inter_pairs;

    prop_POS_sig_STNGPe_pairs = sum(pos_sigSTNGPepairs) / n_inter_pairs
    prop_NEG_sig_STNGPe_pairs = sum(neg_sigSTNGPepairs) / n_inter_pairs


   
    if do_total_xcorr_stats
        [mSTNGPe,sSTNGPe,spSTNGPe,stSTNGPe] = xcorr_stats(total_STNGPe_ccf,total_STNGPe_pairs);
        if do_display
            figure
            bar(xSTNGPe,total_STNGPe_ccf,1)
            axis([-0.1 0.1 min(total_STNGPe_ccf)  max(total_STNGPe_ccf)])
            line([-0.1 0.1],[spSTNGPe spSTNGPe],'Color',[0 0 0]);
        end
    end
end

%%% do auto-correlograms %%%%
if do_acorr
    STN_acf = cell(n_STN,1);
    GPe_acf = cell(n_GPe,1);
    for i = 1:n_STN
        fprintf(1, 'STN autocorrelation loop count %d\n', i);
        [STN_acf{i},xSTNacf,f1,f2,n_pairs] = LIF_xcorr(STN_times{i},STN_times{i}, binsize, [0 time_seconds], maxlag);    
    end
    
    % [mSTN,sSTN,spSTN,stSTN] = xcorr_stats(STN_ccf{i},n_pairs,p);
    for i = 1:n_GPe
        fprintf(1, 'GPe autocorrelation loop count %d\n', i);
        [GPe_acf{i},xGPeacf,f1,f2,n_pairs] = LIF_xcorr(GPe_times{i},GPe_times{i}, binsize, [0 time_seconds], maxlag);    
    end
    
    % [mGPe,sGPe,spGPe,stGPe] = xcorr_stats(GPe_acf{i},n_pairs,p);    
end

%%%%%%%%%%%%%%% coherence analysis %%%%%%%%%%%%%
%%%% currently GPe and STN-GPe only %%%%%%%%%%%%%%%%%%%%%%
if do_coherence
   
    %%% do GPe coherence %%%
    GPe_cohere_res = struct('freqs', [], 'cohere', [], 'errs', [], 'peak_freq', [], 'peak_cohere', [], 'meanC', []);
        
    disp('doing GPe coherence');
    for i = 1:n_GPe_pairs
        fprintf(1, 'GPe coherence loop count %d\n', i);
        
        %% GPe %%% 
        data1 = GPe_data(GPe_pairs(i,1)).times;
        data2 = GPe_data(GPe_pairs(i,2)).times;
        [r1 c] = size(data1);
        [r2 c] = size(data2);
        r = min(r1,r2);
        if r >= Min_spikesC
            [C,phi,f,confC,phierr,Cerr]=coherencypt(data1, data2, tapersC, padC, FsC, freq_rangeC, err_barsC, 0,fscorrC);
            GPe_cohere_res(i).freqs = f;
            GPe_cohere_res(i).cohere = C;
            GPe_cohere_res(i).errs = Cerr;
            GPe_cohere_res(i).meanC = mean(C);

            GPe_err_lo = Cerr(1,:);
            sig_test_GPe_cohere = GPe_err_lo - mean(C);
            [GPe_peak GPe_peak_index] = max(sig_test_GPe_cohere);
            
            if GPe_peak >0
                GPe_cohere_res(i).peak_freq = f(GPe_peak_index);
                GPe_cohere_res(i).peak_cohere = GPe_peak;
            else
                GPe_cohere_res(i).peak_freq = NaN;
                GPe_cohere_res(i).peak_cohere = NaN;
            end
        else
            GPe_cohere_res(loop).freqs = [];
            GPe_cohere_res(loop).cohere = [];
            GPe_cohere_res(loop).errs = [];
            GPe_cohere_res(loop).meanC = [];
            GPe_cohere_res(loop).peak_freq = NaN;
            GPe_cohere_res(loop).peak_cohere = NaN;
        end
    end

    
    
    %%%% do mean GP Coherence %%%%%
    
    GPe_coheres = {GPe_cohere_res.cohere};
    GPe_Cfreqs = {GPe_cohere_res.freqs};
    N_freqs = cellfun('length', GPe_Cfreqs);
    N_non_null_pairsGP = length(N_freqs);
    [max_N_fs, index] = max(N_freqs);
    GPe_mu_cohere = zeros(max_N_fs,1);
    GPe_Cfreqs_base = GPe_Cfreqs{index};
    keyboard
    for loop=1:N_non_null_pairsGP
        if N_freqs(loop) > 0
            if N_freqs(loop) ~= max_N_fs
                Coheres = interp1(GPe_Cfreqs{loop}, GPe_coheres{loop}, GPe_Cfreqs_base);
                Coheres = Coheres';
            else
                Coheres = GPe_coheres{loop};
            end
        else
            Coheres = zeros(max_N_fs, 1);
        end
        GPe_mu_cohere = GPe_mu_cohere + Coheres;
    end
    GPe_mu_cohere = GPe_mu_cohere ./ N_non_null_pairsGP;
    
    sig_GPe_peak_Cfreqs = [GPe_cohere_res.peak_freq];
    GPe_sig_Cfs = sig_GPe_peak_Cfreqs(find(sig_GPe_peak_Cfreqs > 0));
    lowgc = GPe_sig_Cfs(GPe_sig_Cfs < 5);
    GPe_Cs = [GPe_cohere_res.meanC];

    %%%%% plot GPe coherence %%%%%
    if do_display
        
        figure
        subplot(2,1,1)
        hist(GPe_sig_Cfs, 20);
        title('significant maximum frequencies - coherence Gpe');
        xlabel('peak freq');

        subplot(2,1,2)
       
        hist(lowgc, 15);
        title('significant maximum frequencies in range 0-5Hz - coherence GPe');
        xlabel('peak freq');

        figure
	subplot(2,1,1)
        plot(GPe_Cfreqs_base, GPe_mu_cohere);
        title('Mean GPe coherence');
        xlabel('frequency');

        subplot(2,1,2)
        
        hist(GPe_Cs, 40);
	title('histogram of mean coherence values')
        drawnow
    end
    fprintf(1, 'mean GPe coherence %f\n\n', mean(GPe_Cs));

    %%%%% do GPe-STN coherence %%%%%

    STNGPe_cohere_res = struct('freqs', [], 'cohere', [], 'errs', [], 'peak_freq', [], 'peak_cohere', [], 'meanC', []);
   
    disp('doing STN-GPe coherence');
    for i = 1:n_inter_pairs
        fprintf(1, 'coherence loop count %d\n', i);
        
        %% STN-GPe %%% 
        data1 = GPe_data(STNGPe_pairs(i,2)).times;
        data2 = STN_data(STNGPe_pairs(i,1)).times;
        [r1 c] = size(data1);
        [r2 c] = size(data2);
        r = min(r1,r2);
        if r >= Min_spikesC
            [C,phi,f,confC,phierr,Cerr]=coherencypt(data1, data2, tapersC, padC, FsC, freq_rangeC, err_barsC, 0,fscorrC);
            STNGPe_cohere_res(i).freqs = f;
            STNGPe_cohere_res(i).cohere = C;
            STNGPe_cohere_res(i).errs = Cerr;
            STNGPe_cohere_res(i).meanC = mean(C);

            STNGPe_err_lo = Cerr(1,:);
            sig_test_STNGPe_cohere = STNGPe_err_lo - mean(C);
            [STNGPe_peak STNGPe_peak_index] = max(sig_test_STNGPe_cohere);
            
            if STNGPe_peak >0
                STNGPe_cohere_res(i).peak_freq = f(STNGPe_peak_index);
                STNGPe_cohere_res(i).peak_cohere = STNGPe_peak;
            else
                STNGPe_cohere_res(i).peak_freq = NaN;
                STNGPe_cohere_res(i).peak_cohere = NaN;
            end
        else
            STNGPe_cohere_res(loop).freqs = [];
            STNGPe_cohere_res(loop).cohere = [];
            STNGPe_cohere_res(loop).errs = [];
            STNGPe_cohere_res(loop).meanC = [];
            STNGPe_cohere_res(loop).peak_freq = NaN;
            STNGPe_cohere_res(loop).peak_cohere = NaN;
        end
    end

 %%%% do mean STN-GP Coherence %%%%%
    
    STNGPe_coheres = {STNGPe_cohere_res.cohere};
    STNGPe_Cfreqs = {STNGPe_cohere_res.freqs};
    N_freqs = cellfun('length', STNGPe_Cfreqs);
    N_non_null_pairsSTNGP = length(N_freqs);
    [max_N_fs, index] = max(N_freqs);
    STNGPe_mu_cohere = zeros(max_N_fs,1);
    STNGPe_Cfreqs_base = STNGPe_Cfreqs{index};
    for loop=1:N_non_null_pairsSTNGP
        if N_freqs(loop) > 0
            if N_freqs(loop) ~= max_N_fs
                Coheres = interp1(STNGPe_Cfreqs{loop}, STNGPe_coheres{loop}, STNGPe_Cfreqs_base);
                Coheres = Coheres';
            else
                Coheres = STNGPe_coheres{loop};
            end
        else
            Coheres = zeros(max_N_fs, 1);
        end
        STNGPe_mu_cohere = STNGPe_mu_cohere + Coheres;
    end
    STNGPe_mu_cohere = STNGPe_mu_cohere ./ N_non_null_pairsSTNGP;


%%%% plot STN GPe coherences %%%%%

    sig_STNGPe_peak_Cfreqs = [STNGPe_cohere_res.peak_freq];
    STNGPe_sig_Cfs = sig_STNGPe_peak_Cfreqs(find(sig_STNGPe_peak_Cfreqs > 0));
    lowsgc = STNGPe_sig_Cfs(STNGPe_sig_Cfs < 5);
    STN_GPe_Cs = [STNGPe_cohere_res.meanC];

    if do_display
        figure
        subplot(2,1,1)
        hist(STNGPe_sig_Cfs, 20);
        title('significant maximum  frequencies - coherence STN-Gpe');
        xlabel('peak freq');

        subplot(2,1,2)
        hist(lowsgc, 15);
        title('significant maximum  frequencies in range 0-5Hz - coherence STN-GPe');
        xlabel('peak freq');

        figure
	subplot(2,1,1)
        plot(STNGPe_Cfreqs_base, STNGPe_mu_cohere);
        title('Mean STN-GPe coherence')
        xlabel('frequency');

	subplot(2,1,2)
        
        hist(STN_GPe_Cs, 40);
        title('histogram of mean coherence values')
        drawnow
    end

    fprintf(1, 'mean STN-GPe coherence %f\n\n', mean(STN_GPe_Cs));

end

%%%%%%%%%%%%%%%%%%Spike triggered averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_spike_trig_avg

    ctx_times = cell(neurons_per_nucleus,1);
    ctx_IFRbins = cell(neurons_per_nucleus,1);
    ctx_mean = [];
    
    % generate pseudo-EEG using mean of smoothed single-neuron
    % firing rate estimates
    fprintf(1, 'Computing pseudo cortical EEG \n');
    for loop = 1:neurons_per_nucleus
        ctx_times{loop} = t_in_t(in_n == EXT(loop));
        [ctx_IFRbins{loop},ctx_timepoints] = LIF_firingrate(ctx_times{loop}',smooth_win_size,time_seconds,dt,step_size,window);
        if loop==1 ctx_mean = ctx_IFRbins{loop} ./ neurons_per_nucleus; 
        else
            ctx_mean = ctx_mean + ctx_IFRbins{loop} ./ neurons_per_nucleus;  
        end
    end
   

    smooth_steps = length(ctx_timepoints);
    new_dt = step_size * dt;
    win_in_steps = trig_win_half ./ new_dt;
    
    avg_x_axis = linspace(-win_in_steps,win_in_steps,win_in_steps*2+1) .* new_dt;
    
    GPe_avg_wvfrm = zeros(n_GPe,win_in_steps*2+1);
    STN_avg_wvfrm = zeros(n_STN,win_in_steps*2+1);
     

    for loop1 = 1:n_GPe
        GPe_t_stamps = GPe_times{loop1};
        GPe_num_spikes = length(GPe_t_stamps);
          
        
        for loop2 = 1:GPe_num_spikes
            start_t = round(GPe_t_stamps(loop2) ./ new_dt) - win_in_steps;
            end_t = round(GPe_t_stamps(loop2) ./ new_dt) + win_in_steps;
            if start_t > 0 & end_t <= smooth_steps 
                GPe_avg_wvfrm(loop1,:) = GPe_avg_wvfrm(loop1,:) + ctx_mean(start_t:end_t);
            end
        end
        GPe_avg_wvfrm(loop1,:) = GPe_avg_wvfrm(loop1,:) ./ GPe_num_spikes;
    end
    
    for loop1 = 1:n_STN   
        STN_t_stamps = STN_times{loop1};
        STN_num_spikes = length(STN_t_stamps);

        for loop2 = 1:STN_num_spikes
            start_t = round(STN_t_stamps(loop2) ./ new_dt) - win_in_steps;
            end_t = round(STN_t_stamps(loop2) ./ new_dt) + win_in_steps;
            if start_t > 0 & end_t <= smooth_steps 
                STN_avg_wvfrm(loop1,:) = STN_avg_wvfrm(loop1,:) + ctx_mean(start_t:end_t);
            end
        end
        STN_avg_wvfrm(loop1,:) = STN_avg_wvfrm(loop1,:) ./ STN_num_spikes;
    end

%     figure
%     plot(avg_x_axis,GPe_avg_wvfrm(2,:))
%     axis([-win_in_steps*new_dt win_in_steps*new_dt 0 max(ctx_mean)]);
%     drawnow
%     disp('press to continue')
%     pause;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_modes = struct('means', do_means, 'isi', do_isi, 'spectra', do_spectra,...
        'bursts', do_bursts, 'xcorr_intra', do_xcorr_intra, 'xcorr_inter', do_xcorr_inter,...
        'coherence',do_coherence,'spike_trig_avg',do_spike_trig_avg);

an_fname = [pathroot an_fname];

if save_full_analysis
    save(an_fname, 'analysis_modes','t_in_t','ISIbins','limits',...
        'freq_range', 'tapers', 'pad', 'Fs', 'err_bars', 'fscorr', 'binsize','maxlag','p','n_STN','n_GPe','n_GPi');
    if do_means
        save(an_fname, 'mean_STN','mean_GPe','mean_GPi','STN_Hz','GPe_Hz','GPi_Hz',...
            'sem_STN', 'std_STN', 'sem_GPe', 'std_GPe', 'sem_GPi', 'std_GPi', 'CV_STN', 'CV_GPe', 'CV_GPi', '-append');
    end
    if do_isi
        save(an_fname, 'STN_rates','STN_ISIhist','STN_x_isi','STN_x_hist',...
            'GPe_rates','GPe_ISIhist','GPe_x_isi','GPe_x_hist', '-append');
    end
    if do_spectra
        save(an_fname, 'STN_spect_res', 'GPe_spect_res',...
             'STN_sig_fs', 'GPe_sig_fs','-append');
    end
    if do_bursts
        save(an_fname, 'N_STN_bursts', 'N_GPe_bursts',...
            'STN_bursts','STN_burst_isis', 'GPe_bursts','GPe_burst_isis', '-append');
    end
    if do_xcorr_intra
        save(an_fname, 'n_intra_pairs', 'STN_pairs','GPe_pairs', ...
            'STN_ccf','xSTN','mean_STN_ccf','GPe_ccf','xGPe','mean_GPe_ccf',...
            'prop_POS_sig_STN_pairs','prop_NEG_sig_STN_pairs','prop_POS_sig_GPe_pairs','prop_NEG_sig_GPe_pairs',...
            'STN_peak_and_trough','GPe_peak_and_trough', '-append');
    end
    if do_xcorr_inter
        save(an_fname,'n_inter_pairs','STNGPe_pairs',...
            'STNGPe_ccf','xSTNGPe','mean_STNGPe_ccf',... 
            'prop_POS_sig_STNGPe_pairs','prop_NEG_sig_STNGPe_pairs','STNGPe_peak_and_trough','-append');
    end
    if do_acorr 
        save(an_fname,'STN_acf','xSTNacf','GPe_acf','xGPeacf','-append');
    end
    if do_coherence
        save(an_fname, 'n_intra_pairs','GPe_pairs', 'GPe_cohere_res', 'GPe_Cfreqs_base', 'GPe_mu_cohere',...
	     'n_inter_pairs','STNGPe_pairs', 'STNGPe_cohere_res', 'STNGPe_Cfreqs_base', 'STNGPe_mu_cohere', '-append');
    end
   if do_spike_trig_avg
        save(an_fname, 'ctx_mean', 'GPe_avg_wvfrm','STN_avg_wvfrm', '-append');
    end
end
if save_summary_analysis
    
    save(an_fname, 'neurons_per_nucleus', 'maxlag','analysis_modes');

    if do_means
        save(an_fname, 'mean_STN','mean_GPe','mean_GPi','STN_Hz', 'sem_STN', 'std_STN',...
            'GPe_Hz', 'sem_GPe', 'std_GPe', 'GPi_Hz', 'sem_GPi', 'std_GPi', 'CV_STN', 'CV_GPe', 'CV_GPi', '-append');
    end
    if do_isi
        save(an_fname, 'STN_times', 'GPe_times', 'rs_STN', 'rs_GPe', 'rs_in', 'STN_rates', 'STN_x_isi',...
	     'GPe_rates', 'GPe_x_isi', '-append');
    end
    if do_spectra
        save(an_fname, 'STN_freqs_base', 'STN_mu_power', 'GPe_freqs_base', 'GPe_mu_power',...
            'STN_sig_fs', 'GPe_sig_fs', 'lows', 'lowg', '-append');
    end
    if do_bursts
        save(an_fname, 'N_STN_bursts', 'N_GPe_bursts', '-append');
    end
    if do_xcorr_intra 
        save(an_fname,'xSTN', 'total_STN_ccf', 'xGPe', 'total_GPe_ccf', 'n_intra_pairs',...
            'prop_POS_sig_STN_pairs','prop_NEG_sig_STN_pairs','prop_POS_sig_GPe_pairs','prop_NEG_sig_GPe_pairs',...
            'STN_peak_and_trough','GPe_peak_and_trough','-append');
    end
    if do_xcorr_inter
        save(an_fname,'n_inter_pairs', 'total_STNGPe_ccf','xSTNGPe',... 
            'prop_POS_sig_STNGPe_pairs','prop_NEG_sig_STNGPe_pairs','STNGPe_peak_and_trough', '-append');
    end
    if do_acorr 
        save(an_fname,'STN_acf','xSTNacf','GPe_acf','xGPeacf','-append');
    end
    if do_coherence
        save(an_fname, 'n_intra_pairs', 'N_non_null_pairsGP','N_non_null_pairsSTNGP',  'GPe_Cfreqs_base', 'GPe_mu_cohere', 'GPe_Cs', 'STN_GPe_Cs',...
	      'STNGPe_Cfreqs_base', 'STNGPe_mu_cohere', 'GPe_sig_Cfs', 'STNGPe_sig_Cfs',...
	      'lowgc', 'lowsgc', '-append');
    end
    if do_spike_trig_avg
        save(an_fname, 'ctx_mean', 'GPe_avg_wvfrm', 'STN_avg_wvfrm', '-append');
    end
   
end
