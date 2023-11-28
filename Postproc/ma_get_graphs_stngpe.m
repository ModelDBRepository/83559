%%%% script to display graphs from models-as-animals data (LFO single cell fitting)

clear all

%%
exp_no = 18;
freq_display_limit = 40;

%% paths
batch_path = '../ResultsArchive/LFO-urethane/CompleteModel/ConditionA/';
%batch_path = '../ResultsArchive/LFO-urethane/CompleteModel/ConditionB/';
%batch_path = '../ResultsArchive/LFO-urethane/CompleteModel/ConditionC/';
%batch_path = '../ResultsArchive/LFO-urethane/CompleteModel/ConditionD/';

% Collaterals, no STN 3rd DA: batch of 50 - Normal model
batch = 'LFO_a_20060407T004326_batch.mat';
%batch = 'LFO_b_20060407T011449_batch.mat';
%batch = 'LFO_c_20060407T003124_batch.mat'; 
%batch = 'LFO_d_20060407T000419_batch.mat';


   % parameters for analyses done here
    ISIbins = 100; %No. of bins in ISI hist.
    limits = [0.005 0.6];  % ISI hist limits in seconds (outside these bounds get put into end bins)
    trig_win_half = 1;
    acf_binsize = 0.01;   % in seconds
    acf_maxlag = 2;     % window size for auto-correlograms

% parameters for periodogram
freq_range = [0.1 10];    % min and max freqs in periodogram
tapers = [3 5];                         % recommended in Chronuz; 5 tapers with bandwidth*time = 3
pad = 2;                                % Padding the frequency sampling by another 2 powers of 2;
err_bars = [1 0.01];                    % method 1 (theoretical) and p_sig = 0.01
fscorr = 0;                             %don't use finite size corrections

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
load([batch_path batch]);

[n_batches c] = size(batch_analysis_list);

%%%%%%%%%% do Condition A (control) %%%%%%%%%%%%%%
extracted_list = batch_analysis_list{exp_no,3};
analysis_list= batch_analysis_list{exp_no,2};

n_animals = length(analysis_list);

cl_fig
for loop1 = 1:n_animals

     load([batch_path extracted_list{loop1}])
     load([batch_path analysis_list{loop1}])

    fprintf(1, 'mean STN rate %.2f:  stderr of mean STN rate %.2f: std dev of STN rate %.2f\n', STN_Hz, sem_STN, std_STN);
    fprintf(1, 'mean GPe rate %.2f:  stderr of mean GPe rate %.2f: std dev of STN rate %.2f\n', GPe_Hz, sem_GPe, std_GPe);
    fprintf(1, 'mean GPi rate %.2f:  stderr of mean GPi rate %.2f: std dev of STN rate %.2f\n', GPi_Hz, sem_GPi, std_GPi);

    % do cortical EEG power spectrum
    new_dt = dt * 10;
    time_stamps = new_dt:new_dt:time_seconds;
    [ctx_result,ctx_power,ctx_fs] = scargle_analysis(ctx_mean,time_stamps,new_dt,time_seconds,256,[0.1 5]);

    %%%%%%%%%% show data for STN %%%%%%%%%

    STN_acf = cell(n_STN,1);
    xSTNacf = cell(n_STN,1);
    temp_acf = cell(n_STN,1);
    limit = cell(n_STN,1);    

    % do power spectrum of cortical wave
    % Fs = 1 ./ dt;                           % for spectrum  best to use the underlying sampling;
    % [ctx_MT_power, ctx_MT_fs_array, R, ctx_MT_errs] = mtspectrumpt(ctx_mean', tapers, pad, Fs, freq_range, err_bars, 0, fscorr);
 
    for loop = 1:n_STN
        figure
 
        subplot(2,3,1), plot(ctx_mean);
        [STN_rates{loop},STN_ISIhist{loop},STN_x_isi{loop},STN_x_hist] = LIF_ISI_analysis(STN_times{loop},ISIbins,limits);
        subplot(2, 3, 4);
        plot(STN_x_isi{loop}, STN_rates{loop});
       
        subplot(2,3,2), plot(ctx_fs,ctx_power)  
        % run auto-corr, take out centre bin
        [STN_acf{loop},xSTNacf{loop},f1,f2,n_pairs] = LIF_xcorr(STN_times{loop},STN_times{loop}, acf_binsize, [0 time_seconds], acf_maxlag);    
 
	temp_acf{loop} = STN_acf{loop};
        temp_acf{loop}(find(xSTNacf{loop}>=-0.05 & xSTNacf{loop}<= 0.05)) = mean(STN_acf{loop});
       
        
        subplot(2,3,5), plot(xSTNacf{loop},temp_acf{loop})
        
        % plot cell spectra 0-10Hz
	limit{loop} = find(STN_spect_res(loop).freqs <= freq_display_limit);
        subplot(2,3,3), plot(STN_spect_res(loop).freqs(limit{loop}),STN_spect_res(loop).powers(limit{loop}));

        % make time-scale for x-axis (+/- 1 s?)
	trig_x = linspace(-trig_win_half,trig_win_half,length(STN_avg_wvfrm(loop,:)));			
     
        % plot spike-triggered average waveform
        subplot(2,3,6), plot(trig_x,STN_avg_wvfrm(loop,:));

        h = gca;
        txt = ['animal ' num2str(loop1) '- STN, cell ' num2str(loop)];
        title(txt);
        drawnow
        

    end
    tile
    resp = input('save STN cell data for? (enter each number e.g. 41, enter 0 to skip)');
    if resp == 0
       
    else
       resp = num2str(resp);
       for i = 1:length(resp);
          cell_num = str2num(resp(i));
          file_prefix = [batch(1:end-9) '_exp_' num2str(exp_no) '_animal_' num2str(loop1) '_STNcell_' resp(i)];

          % ctx eeg
          ctx_eeg = [time_stamps' ctx_mean'];
          save([file_prefix '_ctx_eeg.txt'],'ctx_eeg','-ascii');

          % generate spikes before saving
          ts = STN_times{cell_num};
          marks = -60 .* ones(length(ts), 1);
          stn_ras = [ts' marks];
          save([file_prefix '_ras.txt'],'stn_ras','-ascii'); 

          % cortex Lomb
          ctx_Lomb = [ctx_fs' ctx_power'];
          save([file_prefix '_ctx_Lomb.txt'],'ctx_Lomb','-ascii');

          % autocorrelogram
          stn_acorr = [xSTNacf{cell_num}' temp_acf{cell_num}'];  
          save([file_prefix '_acorr.txt'],'stn_acorr','-ascii');

          % power spectrum
          stn_power = [STN_spect_res(cell_num).freqs' STN_spect_res(cell_num).powers];
          save([file_prefix '_power.txt'],'stn_power','-ascii');

          % spike-triggered avg waveform
          stn_wvfrm = [trig_x' STN_avg_wvfrm(cell_num,:)'];
          save([file_prefix '_avg_wvfrm.txt'],'stn_wvfrm','-ascii');
       end
    end

    cl_fig

    % and the same for the GPe.....
 
    GPe_acf = cell(n_GPe,1);
    xGPeacf = cell(n_GPe,1);
    temp_acf = cell(n_GPe,1);
    limit = cell(n_GPe,1);
  
    for loop = 1:n_GPe
        figure
 
        subplot(2,3,1), plot(ctx_mean);
        [GPe_rates{loop},GPe_ISIhist{loop},GPe_x_isi{loop},GPe_x_hist] = LIF_ISI_analysis(GPe_times{loop},ISIbins,limits);
        subplot(2, 3, 4);
        plot(GPe_x_isi{loop}, GPe_rates{loop});
       
        subplot(2,3,2), plot(ctx_fs,ctx_power)  
        % run auto-corr, take out centre bin
        [GPe_acf{loop},xGPeacf{loop},f1,f2,n_pairs] = LIF_xcorr(GPe_times{loop},GPe_times{loop}, acf_binsize, [0 time_seconds], acf_maxlag);    
 
	temp_acf{loop} = GPe_acf{loop};
        temp_acf{loop}(find(xGPeacf{loop}>=-0.05 & xGPeacf{loop} <= 0.05)) = mean(GPe_acf{loop});
        subplot(2,3,5), plot(xGPeacf{loop},temp_acf{loop})
        
        % plot cell spectra 0-10Hz
	limit{loop} = find(GPe_spect_res(loop).freqs <= freq_display_limit);
        subplot(2,3,3), plot(GPe_spect_res(loop).freqs(limit{loop}),GPe_spect_res(loop).powers(limit{loop}));

        % make time-scale for x-axis (+/- 1 s?)
	trig_x = linspace(-trig_win_half,trig_win_half,length(GPe_avg_wvfrm(loop,:)));			
     
        % plot spike-triggered average waveform
        subplot(2,3,6), plot(trig_x,GPe_avg_wvfrm(loop,:));

        h = gca;
        txt = ['animal ' num2str(loop1) '- GPe, cell ' num2str(loop)];
        title(txt);
        drawnow
        

    end
    tile
    resp = input('save GPe cell data for? (enter each number e.g. 41, enter 0 to skip)');
    if resp == 0
       
    else
       resp = num2str(resp);
       for i = 1:length(resp);
          cell_num = str2num(resp(i));
          file_prefix = [batch(1:end-9) '_exp_' num2str(exp_no) '_animal_' num2str(loop1) '_GPecell_' resp(i)];

          % ctx eeg
          ctx_eeg = [time_stamps' ctx_mean'];
          save([file_prefix '_ctx_eeg.txt'],'ctx_eeg','-ascii');

          % generate spikes before saving
          ts = GPe_times{cell_num};
          marks = -60 .* ones(length(ts), 1);
          stn_ras = [ts' marks];
          save([file_prefix '_ras.txt'],'stn_ras','-ascii'); 

          % cortex Lomb
          ctx_Lomb = [ctx_fs' ctx_power'];
          save([file_prefix '_ctx_Lomb.txt'],'ctx_Lomb','-ascii');

          % autocorrelogram
          stn_acorr = [xGPeacf{cell_num}' temp_acf{cell_num}'];  
          save([file_prefix '_acorr.txt'],'stn_acorr','-ascii');

          % power spectrum
          stn_power = [GPe_spect_res(cell_num).freqs' GPe_spect_res(cell_num).powers];
          save([file_prefix '_power.txt'],'stn_power','-ascii');

          % spike-triggered avg waveform
          stn_wvfrm = [trig_x' GPe_avg_wvfrm(cell_num,:)'];
          save([file_prefix '_avg_wvfrm.txt'],'stn_wvfrm','-ascii');
       end
    end

    cl_fig
 end
  
 
