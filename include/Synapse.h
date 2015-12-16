//
// Created by Matthew Guay on 9/17/15.
//
#include <deque>

#pragma once

namespace kcnet {
    class CholinergicSynapse {
    public:
        double g, E, Oc, tnow, t0, dt, lambda;
        std::deque<double> ts;
        double weight[4];
        double gOcs[4] = {0., 0., 0., 0.};
        bool has_t0;

        CholinergicSynapse(double g_, double dt_, double t0_, double lambda_=0, bool has_t0_=true);
        ~CholinergicSynapse();

        void update(double t_);

        void generate_random_spikes(double lambda_, unsigned int n_=2);

        double T(double t_);

    };
}
