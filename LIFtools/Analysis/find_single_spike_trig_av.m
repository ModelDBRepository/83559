function [avg_x_axis, sktav] = find_single_spike_trig_av(t_stamps, ctx_mean)

trig_win_half = 1;  % in seconds (should have been read in from analysis but kludged here)
step_size = 10;     % window step in time-steps (kludged - likewise)
dt = 0.0001;        % raw time step (kludge again)

new_dt = step_size * dt;
win_in_steps = trig_win_half ./ new_dt;

num_spikes = length(t_stamps);
smooth_steps = length(ctx_mean);

avg_x_axis = linspace(-win_in_steps,win_in_steps,win_in_steps*2+1) .* new_dt;
sktav = zeros(1,win_in_steps*2+1);
for loop2 = 1:num_spikes
    start_t = round(t_stamps(loop2) ./ new_dt) - win_in_steps;
    end_t = round(t_stamps(loop2) ./ new_dt) + win_in_steps;
    if start_t > 0 & end_t <= smooth_steps
        sktav = sktav + ctx_mean(start_t:end_t);
    end
end
sktav = sktav ./ num_spikes;
