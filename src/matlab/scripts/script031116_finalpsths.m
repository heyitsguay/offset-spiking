%% script031116_finalpsths.m
% Create PSTHs for your mothafuckin' thesiiiisss

%% Load the trials
tr = load_trials('test031116_100');

%% Make PSTHs.
% 20 ms
create_psths(tr, 20);

% 50 ms
create_psths(tr, 50);

% 100 ms
create_psths(tr, 100);

% 130 ms
create_psths(tr, 130);

% 160 ms
create_psths(tr, 160);