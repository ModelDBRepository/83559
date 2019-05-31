% script to summarise output of GHS model i.e. from SNr/GPi nucleus
% Mark Humphries 1/12/2004

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



cl_fig

model = 7;
sim_num = 50;

%file = '../ResultsArchive/Selection/Normal/selection_20060407T015909_batch.mat';
%file = '../ResultsArchive/Selection/LowDA/selection_lowDA_20060407T020953_batch.mat';
file = '../ResultsArchive/Selection/HighDA/selection_highDA_20060407T025808_batch.mat';

load(file) 

% load result file from first batch of simulations
results_file = batch_analysis_list{3}{model}{sim_num};

load(results_file)

t_out_t = double(out_t) .* dt;
t_in_t = double(in_t) .* dt;
time_steps = time_seconds / dt; 

clear link_n link_w

% parameters for firing-rate functions
step_size = 10; % time-steps per windowed value
win_size = 0.05; % seconds (normally 0.05)

%%%%%%%%%%%%%% show raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trace data display
if any(trace_n == GPi)
    trace_times = 1:1:time_steps;
    figure
    plot(trace_times,trace_vals(:,1))
    figure
    plot(trace_times,trace_vals(:,2))
    figure
    plot(trace_times,trace_vals(:,3))
end

% initialise parameters for analysis
GPi_events = cell(neurons_per_nucleus,1);
GPi_times = cell(neurons_per_nucleus,1);
GPi_IFRbins = cell(neurons_per_nucleus,1);


% mean outputs
% mean_GPi = zeros(neurons_per_nucleus,1); 
% 
% for loop = 1:neurons_per_nucleus
%     mean_GPi(loop) = sum(out_n==GPi(loop)) / time_seconds;
% end
% GPi_Hz = mean(mean_GPi);

% raster display
% raster_plot(out_n,out_t);

%%%%% structured_raster(in_n, in_t, out_n, out_t, dt, neurons_per_nucleus, SD1, SD2, STN, GPe, GPi, time_seconds);

% smooth outputs
for loop = 1:neurons_per_nucleus
    GPi_times{loop} = t_out_t(out_n == GPi(loop));
    [GPi_IFRbins{loop},IFR_timepoints] = LIF_firingrate(GPi_times{loop}',win_size,time_seconds,dt,step_size,'alpha');
end


%%%%%%%% do isi histograms %%%%%%%%%%%%%
% parameters for ISI data
ISIbins = 35; %No. of bins in ISI hist.
limits = [0.01 0.6];  % ISI hist limits in seconds (outside these bounds get put into end bins)
No_t_segs = length(switches);

No_chans = 3;

GPi_rates = cell(length(GPi),No_t_segs);            % instantaneous rates for *all* events (not just bursts)
GPi_x_isi = cell(length(GPi),No_t_segs);            % time stamps for rates above (almost just timestamps )
GPi_ISIhist = cell(length(GPi),No_t_segs);          % ISI histograms
GPi_mean_rates = cell(No_chans, No_t_segs);

% GPi_cum_hist = cell(No_chans, No_t_segs);
% for i = 1:No_chans
%     for j = 1:No_t_segs
%         GPi_cum_hist{i,j} = 0;
%     end
% end

max_count = 0;
max_rate = 0;

for loop = 1:neurons_per_nucleus
    t_start = 0;
    for j = 1:No_t_segs
        t_end = switches(j) .* dt;
        duration_seg = t_end - t_start;
        t_times = GPi_times{loop};
        bool_times = t_times > t_start & t_times < t_end;
        times = t_times(bool_times);
        [GPi_rates{loop, j},GPi_ISIhist{loop, j},GPi_x_isi{loop, j},GPi_x_hist] = LIF_ISI_analysis(times,ISIbins,limits);
        chan = floor((loop - 1) ./ neurons_per_channel) + 1;
        if isempty(GPi_rates{loop, j})
            mu_rate = 0;
        else
            mu_rate = mean(GPi_rates{loop, j});
        end
        
        % check for maximum rate
        if mu_rate > max_rate max_rate = mu_rate; end

        GPi_mean_rates{chan, j} = [GPi_mean_rates{chan,j} mu_rate];
        % GPi_cum_hist{chan,j} = GPi_cum_hist{chan,j} + mean(GPi_rates{loop, j});
        t_start = t_end;
    end
end

GPi_mean_dist = cell(No_chans, No_t_segs);          % histograms of GPi output 

No_bins = 20;
bin_edges = linspace(0,ceil(max_rate),No_bins);

for loop1 = 1:No_chans
     for loop2 = 1:No_t_segs
        dist_hist  = histc(GPi_mean_rates{loop1,loop2},bin_edges);
        % find greatest value bin
        if max(dist_hist) > max_count max_count = max(dist_hist); end
        GPi_mean_dist{loop1,loop2} = dist_hist;
     end
end

f_hist = figure
subplot(2,3,1)
h = bar(bin_edges,GPi_mean_dist{1,1}, 'histc');
axis([0 max(bin_edges) 0 max_count]);
p = get(h,'Parent');
set(p,'FontSize',14);

subplot(2,3,2)
h = bar(bin_edges,GPi_mean_dist{1,2}, 'histc');
axis([0 max(bin_edges) 0 max_count]);
p = get(h,'Parent');
set(p,'FontSize',14);

subplot(2,3,3)
h = bar(bin_edges,GPi_mean_dist{1,3}, 'histc');
axis([0 max(bin_edges) 0 max_count]);
p = get(h,'Parent');
set(p,'FontSize',14);

subplot(2,3,4)
h = bar(bin_edges,GPi_mean_dist{2,1}, 'histc');
axis([0 max(bin_edges) 0 max_count]);
p = get(h,'Parent');
set(p,'FontSize',14);

subplot(2,3,5)
h = bar(bin_edges,GPi_mean_dist{2,2}, 'histc');
axis([0 max(bin_edges) 0 max_count]);
p = get(h,'Parent');
set(p,'FontSize',14);

subplot(2,3,6)
h = bar(bin_edges,GPi_mean_dist{2,3}, 'histc');
axis([0 max(bin_edges) 0 max_count]);
p = get(h,'Parent');
set(p,'FontSize',14);


% summarise channel outputs
GPi_ch1 = GPi_IFRbins{1};
GPi_ch2 = GPi_IFRbins{neurons_per_channel+1};
GPi_ch3 = GPi_IFRbins{neurons_per_channel*2+1};

% plot means
for loop1 = 2:neurons_per_channel
    GPi_ch1 = GPi_ch1 + GPi_IFRbins{loop1}; 
    GPi_ch2 = GPi_ch2 + GPi_IFRbins{loop1+neurons_per_channel}; 
    GPi_ch3 = GPi_ch3 + GPi_IFRbins{loop1+neurons_per_channel*2}; 
end
% compute mean for plotting
GPi_ch1 = GPi_ch1 ./ neurons_per_channel;
GPi_ch2 = GPi_ch2 ./ neurons_per_channel;
GPi_ch3 = GPi_ch3 ./ neurons_per_channel;

%%%%%%%% basic mean firing rate stats

% determine mean outputs (firing rates) for each neuron
mean_STN = zeros(neurons_per_nucleus,1);
mean_SD1 = zeros(neurons_per_nucleus,1);
mean_SD2 = zeros(neurons_per_nucleus,1);

for loop = 1:neurons_per_nucleus
    mean_STN(loop) = sum(out_n==STN(loop)) / time_seconds;
    mean_SD1(loop) = sum(out_n==SD1(loop)) / time_seconds;
    mean_SD2(loop) = sum(out_n==SD2(loop)) / time_seconds;
end

% find mean across neurons in each nucleus
STN_Hz = mean(mean_STN);
SD1_Hz = mean(mean_SD1);
SD2_Hz = mean(mean_SD2);

% find std errors of firing rates in each nucleus

std_STN = std(mean_STN);
std_SD1 = std(mean_SD1);
std_SD2 = std(mean_SD2);

sem_STN = std_STN / sqrt(neurons_per_nucleus);
sem_SD1 = std_SD1 / sqrt(neurons_per_nucleus);
sem_SD2 = std_SD2 / sqrt(neurons_per_nucleus) ;       % computes the standard error of the mean (SEM) as quoted in e.g Urbain et al 2000

fprintf(1, 'mean STN rate %.2f:  stderr of mean STN rate %.2f: std dev of STN rate %.2f\n', STN_Hz, sem_STN, std_STN);
fprintf(1, 'mean SD1 rate %.2f:  stderr of mean SD1 rate %.2f: std dev of STN rate %.2f\n', SD1_Hz, sem_SD1, std_SD1);
fprintf(1, 'mean SD2 rate %.2f:  stderr of mean SD2 rate %.2f: std dev of STN rate %.2f\n', SD2_Hz, sem_SD2, std_SD2);

% also compute SD for population responses!!
f_out = figure;
h = plot(IFR_timepoints,GPi_ch1,'LineWidth',1);
hold on
plot(IFR_timepoints,GPi_ch2,'r','LineWidth',1)
plot(IFR_timepoints,GPi_ch3,'k','LineWidth',1)
hold off
p = get(h,'Parent');
set(p,'FontSize',14);

figure
plot(IFR_timepoints,GPi_IFRbins{1})
hold on
for loop = 2:neurons_per_channel
    plot(IFR_timepoints,GPi_IFRbins{loop}) 
end
title('All GPi channel 1 outputs')
hold off

figure
plot(IFR_timepoints,GPi_IFRbins{neurons_per_channel+1})
hold on
for loop = neurons_per_channel+2:2*neurons_per_channel
    plot(IFR_timepoints,GPi_IFRbins{loop}) 
end
title('All GPi channel 2 outputs')

figure
plot(IFR_timepoints,GPi_IFRbins{2*neurons_per_channel+1})
hold on
for loop = 2*neurons_per_channel+2:3*neurons_per_channel
    plot(IFR_timepoints,GPi_IFRbins{loop}) 
end
title('All GPi channel 3 outputs')


tile

%strFile1 = [file '_hist.png'];
%strFile2 = [file '_outputs.png'];

%print(f_hist,'-dpng','-r600',strFile1);
%print(f_out,'-dpng','-r600',strFile2);

bin_edges = bin_edges';
IFR_timepoints = IFR_timepoints';
GPi_ch1 = GPi_ch1';
GPi_ch2 = GPi_ch2';
GPi_ch3 = GPi_ch3';

save('bin_edges.txt','bin_edges','-ascii')
save('IFR_timepoints.txt','IFR_timepoints','-ascii')
save('GPi_ch1_meanISI.txt','GPi_ch1','-ascii')
save('GPi_ch2_meanISI.txt','GPi_ch2','-ascii')
save('GPi_ch3_meanISI.txt','GPi_ch3','-ascii')

GPi_ch1_t1 = GPi_mean_dist{1,1}' ./ neurons_per_channel .* 100;
GPi_ch1_t2 = GPi_mean_dist{1,2}' ./ neurons_per_channel .* 100;
GPi_ch1_t3 = GPi_mean_dist{1,3}' ./ neurons_per_channel .* 100;
GPi_ch2_t1 = GPi_mean_dist{2,1}' ./ neurons_per_channel .* 100;
GPi_ch2_t2 = GPi_mean_dist{2,2}' ./ neurons_per_channel .*100;
GPi_ch2_t3 = GPi_mean_dist{2,3}'./ neurons_per_channel .*100;

save('GPi_ch1_t1hist.txt','GPi_ch1_t1','-ascii')
save('GPi_ch1_t2hist.txt','GPi_ch1_t2','-ascii')
save('GPi_ch1_t3hist.txt','GPi_ch1_t3','-ascii')
save('GPi_ch2_t1hist.txt','GPi_ch2_t1','-ascii')
save('GPi_ch2_t2hist.txt','GPi_ch2_t2','-ascii')
save('GPi_ch2_t3hist.txt','GPi_ch2_t3','-ascii')
