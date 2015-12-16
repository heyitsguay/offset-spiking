%% script101415_viztrials.m     October 14, 2015
% Creates a visualization of the offset-spiking simulation trial data.

% Matthew Guay
% mguay@math.umd.edu
% Applied Mathematics & Statistics, and Scientific Computation
% Department of Mathematics
% University of Maryland, College Park

%%

% Load trial data.
trialname = 'trial101415_long';
trials = load_trials(trialname);

% Possible synaptic activation windows.
window_lengths = [];
for i = 1:length(trials)
    % Window length for this set of trials.
    w = trials{i}{1}.window;
    
    % Add w to window_lengths if it's not in there already.
    if ~any(abs(w-window_lengths) < 1e-10)
        window_lengths(end+1) = w; %#ok<SAGROW>
    end    
end
% Ascending order.
window_lengths = sort(window_lengths);
%% Turned the cell below into a function.
hf = create_psths(trials, 50);
%%

% Number of trials per parameter configuration.
trials_per = length(trials{1});

% Window length to create histograms for. Returns an error if it's not an
% element of window_lengths.
target_window = window_lengths(1);
if ~any(window_lengths == target_window)
    error('target_window is not a valid window length.');
end
    
% Get all trials with window length target_window.
target_trials = {};
for i = 1:length(trials)
    if trials{i}{1}.window == target_window
        target_trials{end+1} = trials{i}; %#ok<SAGROW>
    end
end
n_target_trials = length(target_trials);

% Plot histograms in order of increasing n_synapses. So, get the sorted
% indices of the target_trials.
target_n_synapses = zeros(n_target_trials, 1);
for i = 1:length(target_trials)
    target_n_synapses(i) = target_trials{i}{1}.n_synapses;
end
% Returns the sorted indices.
[~, idxs_sorted] = sort(target_n_synapses);

% Each plot should be x_plot by y_plot pixels. Stack them into n_cols columns.
x_plot = 170;
y_plot = 80;
n_cols = 2;
plots_per_col = ceil(n_target_trials / n_cols);

% Full plot size.
x_full = (x_plot + 10) * n_cols;
y_full = (y_plot + 5) * plots_per_col;

% Create the output figure.
hf = figure('Position', [1, 1, x_full, y_full], 'Color', 'w');

% For each parameter configuration, we're going to create a histogram of KC
% responses, weighted by the number of trials. In other words, construct a
% binned firing rate estimate for the KC under each test condition.
bin = 3; % (ms).
% Use this for result validation - translates the bin edges forward in time
% to see how results depend on the particulars of the binning.
bin_offset = 0;
% Bin edges. Use the 'unique' command in case the right endpoint gets
% repeated.
if bin_offset == 0
    bin_edges = 0:bin:target_trials{1}{1}.runtime;
else
    bin_edges = bin_offset:bin:target_trials{1}{1}.runtime;
end
% Get bin centers
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
    
    % Trials per parameter configuration
    %trials_per = length(target_trials{1});
    
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
    subplot_idx = (y_idx - 1) * n_cols + x_idx;
    ax_subplot = subplot(plots_per_col, n_cols, subplot_idx);
    hold on
    
    % Plot a dotted line at y=0.5 as a reference.
    plot([0 trial.runtime], [0.5 0.5], 'Color', 0.9 * ones(1,3));
    
    % Plot the activation window.
    area(window_x, window_y, 'EdgeAlpha', 0, 'FaceAlpha', 0.15, 'FaceColor', [0.5, 0.25, 0.25]);
    
    % Plot spikes if there are any.
    if n_spikes > 0
        bar(bin_centers, hist_vals, 1);
    end
    
    % Set axes
    axis([0, trial.runtime, 0, 1]);
    set(ax_subplot, 'YTick', [0, 1]);
    
    % Title
    if trial_idx == 1
        subplottitle = sprintf('%d synapses', trial.n_synapses);
    else
        subplottitle = num2str(trial.n_synapses);
    end
    title(subplottitle, 'fontweight', 'normal', 'fontsize', 10);
end

% Give the overall figure a title.
titl = sprintf('KC PSTH, %5.1f ms input window\n %d trials per plot, %2.1f ms bins', target_window, trials_per, bin);
supertitle(titl);
% set(gcf, 'NextPlot', 'add');
% axes;
% h = title(titl, 'fontweight', 'normal');
% ax = gca;
% ax.Visible = 'off';
% set(h, 'Visible', 'on');
