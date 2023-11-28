function plot_single_spect(neuron, spect_res)
low_lim = 10;
spect = spect_res(neuron);
powers = spect.powers;
freqs = spect.freqs;

low_freqs = freqs(freqs<low_lim);
no_freqs = length(low_freqs);
low_powers = powers(1:no_freqs);
plot(low_freqs, low_powers);
