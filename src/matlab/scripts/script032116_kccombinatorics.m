%% KC sparsity combinatorics script.

%% Estimate single KC spike probability.

% Number of PN neurons.
M = 830;

% Number of PN->KC connections, the sources of which arechosen  uniformly 
% at random from the PN population.
K = 415;

% Number of simultaneously-active PN connections required to trigger a
% spike in a KC.
kcThreshold = 100;

% A range of active neuron counts to use.
Ns = 130:170;

% The resulting activation probability for a single KC, one value per
% active neuron count.
Ps = zeros(length(Ns));

% Compute KC activation probabilities.
for i = 1:length(Ns)
    % Active neuron count.
    N = Ns(i);
    
    % Activation probability p computed as a sum over the active incident
    % neuron count k.
    p = 0;
    for k = 100:N
        p = p + nchoosek(N, k) * nchoosek(M-N, K-k) / nchoosek(M, K);
    end
    Ps(i) = p;
end

% Use a log plot to visualize the relationship between N and p.
plot(Ns, log10(Ps))
%% Estimate KC population spike count distribution.
K = 50000; % number of KCs
%p = 4.5e-6; % spike probability
S = 0:20; % spike counts
P = zeros(length(S), 1);
for s = S
    P(s+1) = nchoosek(K, s) * p^s * (1-p)^(K-s);
end

area(S, P);    

