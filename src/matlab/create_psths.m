function hfig = create_psths(trials_, window_, bin_size_, n_cols_, bin_offset_)
% CREATE_PSTHS returns a handle to a figure. This figure displays the
% PSTH of a simulated Kenyon cell in response to repeated trials of
% stimulation from varying numbers of simulated PN->KC synapses. The
% synapses fire uniformly at random within a window of length window_ ms.
% The input trials_ may be either the output of load_trials.m (a simple
% trial cell array), or a cell array containing several such outputs (a
% composite trial cell array).
%
% Filename: create_psths.m
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
% hfig = CREATE_PSTHS(trials_, window_) returns a handle displaying KC
% PSTHs, one for each value of n_synapses tested in the given dataset with
% window length window_. Each PSTH is an average over a number of runs
% determined by the trial data generation procedure - retrieved below as
% the variable 'trials_per'. Uses a default bin size of bin_size_=3 with x
% offset bin_edges_=0, and displays the resulting subplots across n_cols_=2
% columns.
%
% hfig = CREATE_PSTHS(trials_, window_, bin_size_) operates the same as
% previously, using a user-specified bin size during PSTH construction.
%
% hfig = CREATE_PSTHS(trials_, window_, bin_size_, n_cols_) operates the
% same as previously, using a user-specified number of columns to display
% the resulting figure subplots.
%
% hfig = CREATE_PSTHS(trials_, window_, bin_size_, n_cols_, bin_edges_)
% operates the same as previously, using a user-specified x offset for the
% PSTH bin edges, to validate results against binning effects.

% Default parameter values.
if nargin < 3
    bin_size_ = 2;
end
if nargin < 4
    n_cols_ = 2;
end
if nargin < 5
    bin_offset_ = 0;
end

% If true, trials_ is a composite trial cell array, else it's a simple
% trial cell array.
composite_trials = iscell(trials_{1}{1});

% Store the lengths of the synaptic activation windows used in the trials
% stored in trials_.
window_lengths = [];
% If trials_ is a composite trial cell array...
if composite_trials
    for i = 1:length(trials_)
        L = length(trials_{i});
        for j = 1:L
            % Window length for this set of trials.
            w = trials_{i}{j}{1}.window;
            
            % Add w to window_lengths if it's not in there already.
            if ~any(abs(w - window_lengths) < 1e-10)
                window_lengths(end + 1) = w; %#ok<AGROW>
            end
        end
    end
% If trials_ is a simple trial cell array...
else
for i = 1:length(trials_)
    % Window length for this set of trials.
    w = trials_{i}{1}.window;
    
    % Add w to window_lengths if it's not in there already.
    if ~any(abs(w-window_lengths) < 1e-10)
        window_lengths(end+1) = w; %#ok<AGROW>
    end
end

% Given an input activation window length of window_, use the closest value
% in window_lengths. Display a warning if window_ isn't any of the values
% in window_lengths.
[~, target_window_idx] = min(abs(window_lengths - window_));
target_window = window_lengths(target_window_idx);
if target_window ~= window_
    warn = sprintf('Warning: input window length not found, using %5.0f instead.\n', target_window);
    disp(warn);
end

% if composite_trials
%     % Number of trials per parameter configuration. Currently, assume it's
%     % constant across the different simple trials cell arrays within the
%     % composite trial cell array.
%     trials_per = length(trials_{1}{1});
%     
%     % Get all trials with window length target_window.
%     target_trials = {};
%     for i = 1:length(trials_)
%         for j = 1:length(trials_{i})
%             if trials_{i}{j}{1
% else
% Number of trials per parameter configuration.
trials_per = length(trials_{1});

% Get all trials with window length target_window.
target_trials = {};
for i = 1:length(trials_)
    if trials_{i}{1}.window == target_window
        target_trials{end+1} = trials_{i}; %#ok<AGROW>
    end
end
n_target_trials = length(target_trials);
% end

% Plot histograms in order of increasing n_synapses. So, get the sorted
% indices of the target_trials.
% target_n_synapses = zeros(n_target_trials, 1);
% for i = 1:length(target_trials)
%     target_n_synapses(i) = target_trials{i}{1}.n_synapses;
% end
% % Returns the sorted indices.
% [~, idxs_sorted] = sort(target_n_synapses);

% Each plot should be x_plot by y_plot pixels. Stack them into n_cols columns.
x_plot = 170;
y_plot = 80;
plots_per_col = ceil(n_target_trials / n_cols_);

% Full plot size.
x_full = (x_plot + 10) * n_cols_;
y_full = (y_plot + 4) * plots_per_col;

% Create the output figure.
hfig = figure('Position', [1, 1, x_full, y_full], 'Color', 'w');

% For each parameter configuration, we're going to create a histogram of KC
% responses, weighted by the number of trials. In other words, construct a
% binned firing rate estimate for the KC under each test condition.

% Construct the bin edges.
bin_edges = bin_offset_:bin_size_:target_trials{1}{1}.runtime;

% Get bin centers.
bin_centers = zeros(length(bin_edges)-1, 1);
for i = 1:length(bin_centers)
    bin_centers(i) = 0.5 * (bin_edges(i+1) + bin_edges(i));
end

% Create each plot.
x_idx = 1;
y_idx = 0;
trial_idx = 0;
while 1
    trial_idx = trial_idx + 1;
    % End once we've run through all the trials.
    if trial_idx > n_target_trials
        break;
    end
    % Increment the y_idx, or set to 1 and increment x_idx if y_idx >
    % plots_per_col.
    y_idx = y_idx + 1;
    if y_idx > plots_per_col
        y_idx = 1;
        x_idx = x_idx + 1;
    end
    
    % Construct histogram.
    spikes = [];
    for i = 1:trials_per
        trial = target_trials{trial_idx}{i};
        spikes = [spikes; trial.spike_times]; %#ok<AGROW>
    end
    n_spikes = length(spikes);
    if n_spikes > 0
        hist_vals = histcounts(spikes, bin_edges);
        % Weight the values to estimate firing rate.
        hist_vals = hist_vals / trials_per;
    end
    
    % Activation window (x,y) pairs needed for plotting.
    window_x = [0, trial.t0, trial.t0, trial.t0 + trial.window, trial.t0 + trial.window, trial.runtime];
    window_y = [0, 0, 1, 1, 0, 0];
    
    % Select the correct subplot to draw on.
    subplot_idx = (y_idx - 1) * n_cols_ + x_idx;
    ax_subplot = subplot(plots_per_col, n_cols_, subplot_idx);
    hold on
    
    % Plot a dotted line at y=0.5 as a reference.
    plot([0 trial.runtime], [0.5 0.5], 'Color', 0.9 * ones(1,3));
    
    % Plot the activation window.
    %ha = area(window_x, window_y, 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'FaceColor', [0.5, 0.25, 0.25]);
    ha = area(window_x, window_y, 'EdgeColor', 'none', 'FaceColor', [0.5, 0.25, 0.25]);
    drawnow; pause(0.01);
    ha.Face.ColorType = 'truecoloralpha';
    ha.Face.ColorData(4) = 255 * 0.15;
    drawnow; pause(0.01);
    
    % Plot spikes if there are any.
    if n_spikes > 0
        bar(bin_centers, hist_vals, 1);
    end
    
    % Set axes.
    axis([0, trial.runtime, 0, 1]);
    set(ax_subplot, 'YTick', [0, 1]);
    set(ax_subplot, 'XTick', [0, trial.runtime]);
    
    % Set title.
    subplottitle = sprintf('%d active synapses', trial.n_active_synapses);
    title(subplottitle, 'fontweight', 'normal', 'fontsize', 10);
    
end

% Give the overall figure a title.
titl = sprintf('KC PSTH, %5.1f ms input window\n %d trials per plot, %2.1f ms bins', target_window, trials_per, bin_size_);
supertitle(titl);

end