%% script121115_noisypsths.m
% Create PSTHs for noisy spiking activity using 20 trials.

%% Load one noisy and one noiseless trial.
tr_noisy = load_trials('test121515_8x_1');
tr = load_trials('test121515_8x_1_0');

%% Make PSTHs for the noisy and noiseless trials.
% 10 ms
create_psths(tr, 10);
create_psths(tr_noisy, 10);

% 20 ms
create_psths(tr, 20);
create_psths(tr_noisy, 20);