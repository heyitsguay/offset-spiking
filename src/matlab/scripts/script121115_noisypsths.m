%% script121115_noisypsths.m
% Create PSTHs for noisy spiking activity using 100 trials.

%% Load the trials in 5 parts.
t1 = load_trials('test121115_8x_1');
t2 = load_trials('test121115_8x_2');
t3 = load_trials('test121115_8x_3');
t4 = load_trials('test121115_8x_4');
t5 = load_trials('test121115_8x_5');

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