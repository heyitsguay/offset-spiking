//
// Created by Matthew Guay on 9/17/15.
//

#ifndef OFFSET_SPIKING_CHOLINERGICSYNAPSE_H
#define OFFSET_SPIKING_CHOLINERGICSYNAPSE_H
#pragma once

#include <deque>

namespace kcnet {
    class CholinergicSynapse {
    public:
        double g, E, Oc, tnow, t0, dt, lambda;
        // Spike times.
        std::deque<double> ts;
        // RK4 update weights.
        double weight[4];
        // RK4 Oc intermediate values, multiplied by g.
        double gOcs[4] = {0., 0., 0., 0.};
        // If true, this synapse simulates an active (as opposed to random)
        // PN spike.
        bool has_t0;

        CholinergicSynapse(double g_, double dt_, double t0_, double lambda_=0, bool has_t0_=true);
        ~CholinergicSynapse();

        void update(double t_);

        void generate_random_spikes(double lambda_, unsigned int n_=2);

        double T(double t_);

    };
}

#endif //OFFSET_SPIKING_CHOLINERGICSYNAPSE_H