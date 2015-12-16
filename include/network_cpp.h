//
// Created by Matthew Guay on 9/18/15.
//
#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>

#include "cell_cpp.h"
#include "synapse_cpp.h"

namespace kcnet {
    struct trial {
        double runtime;
        unsigned int n_synapses, n_active_synapses;
        std::vector<double> synapse_times;
        std::vector<double> spike_times;
        double t0;
        double window;
    };

    class CppNetwork {
    private:
        bool initialized = false;

    public:
        unsigned int n_synapses, n_active_synapses;
        double dt, t0, t1, sigma_noise, syn_weight, t_sim, lambda;
        CppKenyonCell* kc;
        std::vector<CppCholinergicSynapse*> pnkcs;
        std::vector<double> pnkc_ts;
        std::vector<double> ts, Vs, Cas, I_Ls, I_KLs, I_Cas, I_KCas, I_KAs, I_Nas, I_Ks, I_syns, I_noises;

        // Contains results from each run of the CppNetwork.
        std::vector<trial> trials;

        CppNetwork();
        ~CppNetwork();

        void setup(int n_active_synapses_, double dt_, double t0_, double t1_, double sigma_noise_, double syn_weight_=15., double lambda_=0);

        double gsyn_icdf(double i_r);

        void create_synapses(int i_n_synapses, double i_syn_weight=15., bool has_t0_=true);

        void resample_synapses(double i_syn_weight=15.);

        void run(double _runtime, bool _reset=false);

        std::vector<double> count_spikes(double i_threshold=40.);
    };
}