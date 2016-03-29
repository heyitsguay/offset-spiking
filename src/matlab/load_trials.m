function trials = load_trials(trialname_)
% LOAD_TRIALS returns a cell array containing trial structs - data output
% from run_test0 defined in ../src/offset_test.cpp. Trial data is read from
% the save file named trialname_, with save directory assumed to be the
% default ../src/bin/save/ (for now). 
%
% Filename: create_trials.m
% ========
% Created: 10/14/2015
% =======
% Modified: 10/14/2015 "Created"
% ========
% Author: Matthew Guay
% ======  mguay@math.umd.edu
%         Applied Mathematics & Statistics, and Scientific Computation
%         Department of Mathematics
%         University of Maryland, College Park
%
% Usage:
% =====
% trials = LOAD_TRIALS(trialname_) returns cell array 'trials' read from
% the save file ../../data/[trialname_].dat. Each instance of
% run_test0 performs (trials_per) simulations per parameter configuration.
% Therefore, trials is organized into (n_trials / trials_per) cells, each
% containing (trials_per) trial structs. 

savedir = '../../data/';
trialid = fopen([savedir trialname_ '.dat']);

% Skip the first line, it's a comment.
fgetl(trialid);

% Second line has number of trials and trials per parameter configuration.
line2 = fgetl(trialid);
line2split = strsplit(line2);
n_trials = str2num(line2split{1});
trials_per = str2num(line2split{2});

% Number of different test parameter configurations.
n_configs = round(n_trials / trials_per);

% Cell array containing all the trial data.
trials = cell(n_configs, 1);
% Parse the trials
tic
for i = 1:n_configs
    trials{i} = cell(trials_per, 1);
    for j = 1:trials_per
        trials{i}{j} = create_trial(trialid);
    end
end
toc

% Close the file once all data is read in.
fclose(trialid);
end