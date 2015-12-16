//
// Created by Matthew Guay on 9/18/15.
//
#include <cmath>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>

#include "Network.h"

namespace kcnet {
    const int init_n_syns = 200; // Initial reserved size for vector pnkc

    Network::Network(){}
    Network::~Network() {
        if(initialized) {
            delete (kc);
            for (int i = n_synapses - 1; i >= 0; i--) {
                delete (pnkcs[i]);
            }
        }
    }

    void Network::setup(int n_active_synapses_, double dt_, double t0_, double t1_, double sigma_noise_, double syn_weight_, double lambda_) {

        dt = dt_; // Network update time-step

        // Internal simulation clock
        t_sim = 0.;

        // The PN->KC synapses fire at times distributed uniformly at random in the interval [t0, t1].
        t0 = t0_;
        t1 = t1_;

        // KC current noise is normal with standard deviation = sigma_noise.
        sigma_noise = sigma_noise_;

        // Synaptic weights are drawn from an empirical distribution based on (Jortner et al., 2007).
        // The output is scaled by syn_weight to match the observed EPSP distribution.
        syn_weight = syn_weight_;

        // Random synaptic firing rate parameter
        lambda = lambda_;

        // Erase the list of previous trials.
        trials = {};

        // pnkcs if n_synapses_ != n_synapses. // If this Network has been initialized before:
        if(initialized) {
            // Clear all the data storage vectors.
            ts.clear();
            Vs.clear();
            Cas.clear();
            I_Ls.clear();
            I_KLs.clear();
            I_Cas.clear();
            I_KCas.clear();
            I_KAs.clear();
            I_Nas.clear();
            I_Ks.clear();
            I_syns.clear();
            I_noises.clear();
            pnkc_ts.clear();

            // Erase previous synapses, then clear pnkcs.
            if(pnkcs.size() > 0) {
                for (int i = (int)pnkcs.size() - 1; i >= 0; i--) {
                    delete pnkcs[i];
                }
                pnkcs.clear();
            }

            // Reset n_synapses and n_active_synapses.
            n_synapses = 0;
            n_active_synapses = 0;

            create_synapses(n_active_synapses_, syn_weight);

            /*// first check if n_synapses is changing.
            if(n_active_synapses_ != n_active_synapses) {

                // If the new setup has fewer CholinergicSynapses:
                if(n_synapses_ < n_synapses) {
                    // Free allocated memory.
                    for(int i=n_synapses-1; i>=n_synapses_; i--){
                        delete pnkcs[i];
                    }

                    // Erase the vector elements.
                    pnkcs.erase(pnkcs.end() - (n_synapses - n_synapses_), pnkcs.end());

                    // Update n_synapses.
                    n_synapses = n_synapses_;

                // Else, the new setup has more synapses.
                } else {

                    int n_new_synapses = n_synapses_ - n_synapses;
                    create_synapses(n_new_synapses, syn_weight);
                }
            }*/

            // Resample the synapses' activation times and weights.
            resample_synapses(syn_weight);

        } else {
            // This Network has not been initialized yet.

            // Number of initial synapses in the Network.
            n_active_synapses = 0;
            n_synapses = 0;

            // Create the KenyonCell.
            kc = new KenyonCell(dt);

            // pnkcs is a vector of CholinergicSynapses. Reserve memory for init_n_syns initially.
            pnkcs.reserve(init_n_syns);

            // Create in_synapses_ CholinergicSynapses.
            create_synapses(n_active_synapses_, syn_weight);

            // Mark this Network initialized.
            initialized = true;
        }
    }

    double Network::gsyn_icdf(double i_r) {
        // Used Python to build up a discretization of gysn's log10 empirical distribution's ICDF.
        // Interpolate linearly between the given values.

        // Number of samples in the ICDF discretization.
        static const int n_samples = 33;

        // x coordinates.
        static const double xvals[n_samples] = {-4.62320408, -4.58567959, -4.5481551,  -4.51063061, -4.47310612, -4.43558163,
         -4.39805714, -4.36053265, -4.32300816, -4.28548367, -4.24795918, -4.21043469,
         -4.1729102,  -4.13538571, -4.09786122, -4.06033673, -4.02281224, -3.98528776,
         -3.94776327, -3.91023878, -3.87271429, -3.8351898,  -3.79766531, -3.76014082,
         -3.72261633, -3.68509184, -3.64756735, -3.61004286, -3.57251837, -3.53499388,
         -3.49746939, -3.4599449,  -3.42242041};

        // y coordinates
        static const double yvals[n_samples] = {0.0, 0.0007689212727664954, 0.003728083804029162, 0.008908838375890547,
                                        0.016311184988350897, 0.027806263191566136, 0.04797832375316616, 0.07703675247716085,
                                        0.11498154936355161, 0.16082534017207048, 0.21135915862185464, 0.26633616115283515,
                                        0.3257563477650133, 0.3892406282500822, 0.455137407309387, 0.5232543292298755,
                                        0.5935913940115497, 0.6647676774252416, 0.7284533159224077, 0.7832525862571784,
                                        0.8291654884295536, 0.8663842527577259, 0.8965593978810916, 0.9200697668826614,
                                        0.9369153597624345, 0.947466441860511, 0.9565364625981884, 0.9656064833358657,
                                        0.9746765040735428, 0.9836767750717671, 0.9910802095342763, 0.9962635011356, 1.};//0.9992266498757383};

        // Find the two successive elements of yvals which i_r is in between.
        int idx1 = 0;
        while(i_r > yvals[idx1]) {idx1 += 1;}
        int idx0 = idx1 - 1;
        double y0 = yvals[idx0];
        double y1 = yvals[idx1];
        double dy = y1 - y0;

        double d0 = (i_r - y0) / dy;
        double d1 = (y1 - i_r) / dy;

        return d1 * xvals[idx0] + d0 * xvals[idx1];

    }

    void Network::create_synapses(int n_synapses_, double isyn_weight_, bool has_t0_) {
        // Creates in_synapses_ CholinergicSynapses with synaptic distributions scaled by isyn_weight_.
        if(n_synapses_ > 0) {

            // Create the CholinergicSynapses.
            for(int i=0; i<n_synapses_; i++) {
                double log_gsyn, gsyn, t;

                // Generate gsyn value.
                log_gsyn = gsyn_icdf((double)std::rand() / RAND_MAX);
                gsyn = isyn_weight_ * std::pow(10., log_gsyn);

                if(has_t0_) {
                    // Generate spike time value in the window [t0, t1].
                    t = t0 + (t1 - t0) * (double) std::rand() / RAND_MAX;
                } else t = 10000000.;

                // Add a new CholinergicSynapse with the given parameters to pnkcs.
                CholinergicSynapse* new_syn = new CholinergicSynapse(gsyn, dt, t, lambda, has_t0_);
                pnkcs.push_back(new_syn);
                n_synapses += 1;
                if(has_t0_) {
                    pnkc_ts.push_back(t);
                    n_active_synapses += 1;
                }

            }

            std::sort(pnkc_ts.begin(), pnkc_ts.end());
        }
    }

    void Network::resample_synapses(double isyn_weight_) {
        // Resamples the random synapse firing times (same distribution - uniform on [t0, t1]).
        // Also
        if(n_synapses > 0) {
            // Clear pnkc_ts. It will be repopulated below.
            pnkc_ts.clear();
            // Resample each synapse.
            for(int i=0; i<n_synapses; i++) {
                double log_gsyn, gsyn, t;

                // Resample synaptic strength, gsyn.
                log_gsyn = gsyn_icdf((double)std::rand() / RAND_MAX);
                gsyn = isyn_weight_ * std::pow(10., log_gsyn);
                pnkcs[i]->g = gsyn;

                // Resample t0, if pnkcs[i]->has_t0.
                if(pnkcs[i]->has_t0) {
                    t = t0 + (t1-t0) * (double) std::rand() / RAND_MAX;
                    pnkcs[i]->t0 = t;
                    pnkc_ts.push_back(t);
                }

                // Resample random spiking if lambda > 0.
                if(lambda > 0) {
                    pnkcs[i]->generate_random_spikes(lambda, 2);
                }

            }

            // Sort pnkc_ts from smallest to largest.
            std::sort(pnkc_ts.begin(), pnkc_ts.end());
        }
    }

    void Network::run(double runtime_, bool reset_) {
        if(reset_) {
            kc->reset_state();
            resample_synapses(syn_weight);
        }
        // Reset internal clock.
        t_sim = 0.;

        // Synapse reversal potential, currently constant across synapses.
        double syn_E;

        // RK4 update information about synapse g*Oc states.
        std::vector<double> syn_gOcs;

        // Number of simulation steps.
        const unsigned int n_steps = unsigned(int(std::ceil(runtime_ / dt)));

        if(n_synapses > 0) { // If there are any CholinergicSynapses
            // Allocate space in syn_gOcs to store their update data.
            syn_gOcs.resize(4 * n_synapses);
            // Get their reversal potential from the first CholinergicSynapse.
            syn_E = pnkcs[0]->E;
        } else { // No Synapses. Initialize syn_gOcs to hold one element.
            syn_gOcs.push_back(0.);
            syn_E = 0.;
        }

        // Resize vectors of KenyonCell property values at each time step of the simulation.
        ts.resize(n_steps);
        Vs.resize(n_steps);
        Cas.resize(n_steps);
        I_Ls.resize(n_steps);
        I_KLs.resize(n_steps);
        I_Cas.resize(n_steps);
        I_KCas.resize(n_steps);
        I_KAs.resize(n_steps);
        I_Nas.resize(n_steps);
        I_Ks.resize(n_steps);
        I_syns.resize(n_steps);

        // Pre-generate the current noise.
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> current_noise(0., sigma_noise);
        for(int i=0; i<n_steps; i++) {
            I_noises.push_back(current_noise(gen));
        }

        // Update the CholinergicSynapses and KenyonCells
        for(int i=0; i<n_steps; i++) {
            if(n_synapses > 0) {
                // Update CholinergicSynapses.
                for(int j=0; j<n_synapses; j++) {
                    // Run the update.
                    pnkcs[j]->update(t_sim);
                    // Save the CholinergicSynapse's gOc information.
                    for(int k=0; k<4; k++) {
                        int idx = j + k * n_synapses;
                        syn_gOcs[idx] = pnkcs[j]->gOcs[k];
                    }
                }

                // Update the KenyonCell.
                kc->update(n_synapses, syn_gOcs, syn_E, I_noises[i]);

                // Update the KC property value lists.
                ts[i] = t_sim;
                Vs[i] = kc->V;
                Cas[i] = kc->Ca;
                I_Ls[i] = kc->I_L;
                I_KLs[i] = kc->I_KL;
                I_Cas[i] = kc->I_Ca;
                I_KCas[i] = kc->I_KCa;
                I_KAs[i] = kc->I_KA;
                I_Nas[i] = kc->I_Na;
                I_Ks[i] = kc->I_K;
                I_syns[i] = kc->I_syn;

                // Update simulation time.
                t_sim += dt;

            }
        }

        // Create trials here.
        trial new_trial;
        new_trial.runtime = runtime_;
        new_trial.n_synapses = n_synapses;
        new_trial.n_active_synapses = n_active_synapses;
        new_trial.spike_times = count_spikes();
        new_trial.synapse_times = pnkc_ts;
        new_trial.t0 = t0;
        new_trial.window = t1 - t0;
        trials.push_back(new_trial);
    }

    std::vector<double> Network::count_spikes(double i_threshold) {
        bool spiking_now = false;
        double peak_now = -10000.;
        double peak_t = -1.;
        std::vector<double> spike_times;

        for(int i=0; i<Vs.size(); i++) {
            if(Vs[i] > i_threshold && !spiking_now) {
                spiking_now = true;
                peak_now = Vs[i];
                peak_t = i * dt;
            }
            else if(Vs[i] > i_threshold && spiking_now) {
                if(Vs[i] > peak_now) {
                    peak_now = Vs[i];
                    peak_t = i * dt;
                }
            }
            else if(Vs[i] <= i_threshold && spiking_now) {
                spiking_now = false;
                spike_times.push_back(peak_t);
            }
        }

        return spike_times;
    }
}
