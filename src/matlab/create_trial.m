function trial = create_trial(trialid_)
% CREATE_TRIAL returns a struct containing the same fields and data as the
% C++ struct trial defined at the beginning of ../src/network_cpp.h. This
% function assumes trialid_ points to an open file containing trial
% data formatted by the function run_test0 defined in
% ../src/offset_test.cpp. 
%
% Filename: create_trial.m
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
% trial = CREATE_TRIAL(trialid_) returns struct 'trial' read from an open
% file with ID trialid_.

% Initialize the output struct.
trial = struct('runtime', 0, 'n_synapses', 0, 'n_active_synapses', 0, 'synapse_times', 0, 'spike_times', 0, 't0', 0, 'window', 0);

% Each trial is written as 8 lines in the data file. First line is a
% comment, ignore it.
fgetl(trialid_);

% Second line: runtime
line = fgetl(trialid_);
% Value is in ms.
trial.runtime = str2num(line); %#ok<*ST2NM>

% Third line: n_synapses
line = fgetl(trialid_);
trial.n_synapses = str2num(line);

% Fourth line: n_active_synapses
line = fgetl(trialid_);
trial.n_active_synapses = str2num(line);

% Fifth line: synapse_times
line = fgetl(trialid_);
% synapse_times is a vector of n_synapses activation times (ms).
trial.synapse_times = zeros(trial.n_synapses, 1);
split = strsplit(line, ',');
for i = 1:length(split)
    trial.synapse_times(i) = str2num(split{i});
end

% Sixth line: spike_times
% spike_times is either the string 'none' if no spikes were recorded, or a
% vector of spike times in ms.
line = fgetl(trialid_);
% Check to see if any spike times were recorded.
if strcmp(line, 'none')
    trial.spike_times = [];
else
    split = strsplit(line, ',');
    n_spikes = length(split);
    trial.spike_times = zeros(n_spikes, 1);
    for i = 1:n_spikes
        trial.spike_times(i) = str2num(split{i});
    end
end

% Seventh line: t0
line = fgetl(trialid_);
trial.t0 = str2num(line);

% Eighth line: window
line = fgetl(trialid_);
trial.window = str2num(line);

end