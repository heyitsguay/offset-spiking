%% script102815_epsphistogram.m
% Creates a histogram of C++-generated EPSP data.

% Matthew Guay
% mguay@math.umd.edu
% Applied Mathematics & Statistics, and Scientific Computation
% Department of Mathematics
% University of Maryland, College Park

%% Load EPSP data 

% EPSP data file name and location.
filename = 'hist1028_2';
savedir = '../src/bin/save/';

% Open file.
fileid = fopen([savedir filename '.dat']);

% String with all EPSP data. Split at commas.
epsps = strsplit(fgetl(fileid), ',');

% MATLAB array of EPSP data.
epsp_array = zeros(length(epsps), 1);

% Fill the MATLAB array.
for i = 1:length(epsps)
    epsp_array(i) = str2double(epsps{i});
end

% Display a histogram with 7 bins, as in (Jortner et al., 2007). Should be
% similar to the left-hand histogram in Figure 6A of that paper.
hist(epsp_array, 7);
xlabel('$$\log_{10}V$$', 'Interpreter', 'Latex');
titletext = sprintf('EPSP distribution histogram, %d trials', length(epsps));
title(titletext, 'fontweight', 'normal');
set(gcf, 'color', 'w');