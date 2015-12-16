//
// Created by mattguay on 10/11/15.
//
#include "offset_test.h"

#include <omp.h>
#include <fstream>
#include <sys/time.h>
#include <iostream>

#define VERBOSE_SAVES false

namespace kcnet {
    void offset_test(std::string save_name_, const unsigned int n_trials_, const unsigned int n_threads_, const bool use_noise_) {
        struct timeval clock0, clock1;
        gettimeofday(&clock0, NULL);
        std::vector<trial> all_trials;
        const int total_synapses = 400;
        const std::vector<int> Ns = {50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
        const std::vector<double> windows = {10., 20., 30., 40., 50., 60., 100., 115., 130., 145., 160., 175., 190.,
                                             205., 220., 235., 250., 265., 280., 295., 310., 325., 340., 355., 370.,
                                             385., 400.};
        std::vector<std::pair<int, double>> trial_params;
        for (int i = 0; i < Ns.size(); i++) {
            for (int j = 0; j < windows.size(); j++) {
                trial_params.push_back(std::make_pair(Ns[i], windows[j]));
            }
        }

        double lambda = (use_noise_? 0.0025 : 0.);
        double t0 = 10.;
        double dt = 0.03;
        double sigma_noise = 1.5;
        double syn_weight = 15.;

        const unsigned int n_pairs = (const unsigned int) trial_params.size();

        if(n_threads_ == 1) {
            for(int i=0; i<n_pairs; i++) {

                Network net = Network();

                // Get the current pair of trial parameters.
                std::pair<int, double> a_pair = trial_params[i];

                // First parameter: N - number of synapses.
                int N = a_pair.first;

                // Second parameter: window - length of the synapse activation time window.
                double window = a_pair.second;

                // Run for an additional 20 ms to capture spiking at the end of the synapse activation window.
                double runtime = t0 + window + 20.;

                // Set up net.
                net.setup(N, dt, t0, t0 + window, sigma_noise, syn_weight, lambda);
                // Add random-firing synapses if lambda > 0.
                if(lambda > 0) {
                    net.create_synapses(total_synapses - N, syn_weight, false);
                }
                // Set net's KenyonCell to have the desired reset state (at equilibrium).
                //if(i == 0) {
                    net.kc->load_state("eq093015");
                //}

                // Run n_trials_ trials with each parameter pair.
                for (int k = 0; k < n_trials_; k++) {
                    // Run from a previously-computed equilibrium state (i.e. reset the state each time). Random
                    // variables are resampled.
                    net.run(runtime, true);
                }
                all_trials.insert(all_trials.end(), net.trials.begin(), net.trials.end());
            }
//            delete net;
        } else {
            std::vector<trial> *thread_trials = new std::vector<trial>[n_threads_];
            Network *thread_nets = new Network[n_threads_];
//            std::vector<Network> thread_nets;
            for(int i=0; i<n_threads_; i++){
                Network new_net = Network();
                thread_nets[i] = new_net;
                // Run a dummy setup so that we can get the right KenyonCell reset state loaded.
                thread_nets[i].setup(1, dt, t0, t0+1., sigma_noise, syn_weight, lambda);
                thread_nets[i].kc->load_state("eq093015");
            }

            #pragma omp parallel num_threads(n_threads_)
            {
                 #pragma omp for
                 for (int i=0; i<n_pairs; i++) {
                     int thread_idx = omp_get_thread_num();

                     std::pair<int, double> a_pair = trial_params[i];
                     int N = a_pair.first;
                     double window = a_pair.second;
                     double runtime = t0 + window + 20.;

                     Network &net = thread_nets[thread_idx];
                     net.setup(N, dt, t0, t0 + window, sigma_noise, syn_weight, lambda);
                     if(lambda > 0) {
                         net.create_synapses(total_synapses - N, syn_weight, false);
                     }

                     for (int k = 0; k <n_trials_; k++) {
                         net.run(runtime, true);
                     }

                     thread_trials[thread_idx].insert(thread_trials[thread_idx].end(), net.trials.begin(), net.trials.end());
                 }
            }

             for(int i=0; i<n_threads_; i++){
                 all_trials.insert(all_trials.end(), thread_trials[i].begin(), thread_trials[i].end());
             }

//            thread_nets.erase(thread_nets.begin(), thread_nets.end());

            // Delete the created array of vectors of trials.
            delete[] thread_trials;
            delete[] thread_nets;
        }

        // All trials complete at this point.

        // Write trial data to a save file.
        // Format:
        // * Trial [trial #]
        // [runtime]
        // [n_synapses]
        // [n_active_synapses]
        // [synapse_times] (not including random ones)
        // [spike_times]
        // [t0]
        // [window]

        // Open file to write.
        std::ofstream trial_file;
        std::string fname = "../save/" + save_name_ + ".dat";
        trial_file.open(fname.c_str());

        // Quick little header with number of trials and trials per parameter configuration.
        trial_file << "# Header info: (number of trials) (trials per parameter configuration).\n";
        trial_file << all_trials.size() << " " << n_trials_ << "\n";
        // Write trial data.
        for(int i=0; i<all_trials.size(); i++) {
            trial a_trial = all_trials[i];
            // Start with a trial ID.
            trial_file << "* Trial " << i << "\n";
            // Runtime.
            if(VERBOSE_SAVES) trial_file << "Runtime (ms): ";
            trial_file << a_trial.runtime << "\n";
            // Number of synapses.
            if(VERBOSE_SAVES) trial_file << "Number of synapses: ";
            trial_file << a_trial.n_synapses << "\n";
            // Number of actively-firing synapses.
            if(VERBOSE_SAVES) trial_file << "Number of active synapses: ";
            trial_file << a_trial.n_active_synapses << "\n";
            // Synapse activation times.
            if(VERBOSE_SAVES) trial_file << "Synapse activation times (ms): ";
            for(int j=0; j<a_trial.synapse_times.size(); j++) {
                trial_file << a_trial.synapse_times[j];
                if (j < a_trial.synapse_times.size() - 1) {
                    trial_file << ",";
                } else trial_file << "\n";
            }
            // KC spike times.
            if(VERBOSE_SAVES) trial_file << "KC spike times (ms): ";
            if(a_trial.spike_times.size() == 0) {
                // Print a dummy line if no spikes were recorded.
                trial_file << "none\n";
            } else {
                for (int j = 0; j < a_trial.spike_times.size(); j++) {
                    trial_file << a_trial.spike_times[j];
                    if (j < a_trial.spike_times.size() - 1) {
                        trial_file << ",";
                    } else trial_file << "\n";
                }
            }
            // Synapse activation t0.
            if(VERBOSE_SAVES) trial_file << "Synapse activation window begin (ms): ";
            trial_file << a_trial.t0 << "\n";
            // Synapse activation window length.
            if(VERBOSE_SAVES) trial_file << "Synapse activation window length (ms): ";
            trial_file << a_trial.window << "\n";
        }
        trial_file.close();

        gettimeofday(&clock1, NULL);
        double elapsed = ((clock1.tv_sec - clock0.tv_sec) * 1000000u +
                        clock1.tv_usec - clock0.tv_usec) / 1.e6;
        printf("Wrote %d trials to file %s.dat in %f seconds.\n", (int)all_trials.size(), save_name_.c_str(), elapsed);
    }
}
