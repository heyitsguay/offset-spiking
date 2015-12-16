//
// Created by Matthew Guay on 9/17/15.
//
#include "synapse_cpp.h"

#include <vector>
#include <random>
#include <iostream>

namespace kcnet {
    // Quick little Heaviside step function.
    inline double step(double x) {
        return (x < 0)? 0. : 1.;
    }

    // Used in the ODE that updates Oc.
    const double alpha = 0.94;
    const double beta = 0.18;

    CppCholinergicSynapse::CppCholinergicSynapse(double g_, double dt_, double t0_, double lambda_, bool has_t0_) {
        g = g_; // maximum synaptic conductance (μS)
        dt = dt_; // RK4 update time-step (ms)
        t0 = t0_; // Deterministic synaptic activation time.
        lambda = lambda_; // Random activation rate parameter.
        has_t0 = has_t0_; // If true, this synapse participates in coordinates spiking at t0.
        tnow = 0.; // synaptic activation times, relative to simulation start (ms)
        E = 0.; // synaptic reversal potential (mV)
        Oc = 0.; // Proportion of open synaptic channels

        // RK4 steps' time weights.
        weight[0] = 0.;
        weight[1] = 0.5 * dt;
        weight[2] = 0.5 * dt;
        weight[3] = dt;

        // If lambda > 0, generate first two random spike times. Construct ts. Else, just return a deque containing t0.
        generate_random_spikes(lambda);

        // Initialize tnow.
//        tnow = ts.front();
//        ts.pop_front();
    }

    CppCholinergicSynapse::~CppCholinergicSynapse() {}

    void CppCholinergicSynapse::update(double t_) {
        // RK4 update of the synapse's Oc parameter.

        // Update tnow, 0.3 ms after the previous activation.
        if(t_ - tnow > 0.3) {
            // Update tnow if any more ts remain.
            if(ts.size() > 0) {
                tnow = ts.front();
                ts.pop_front();
            }
        }

        double k_Oc[4] = {0., 0., 0., 0.};

        // Rk4 update steps.
        for(int i=0; i<4; i++) {
            double w_i, t_i, Oc_i;
            int iprev = (i == 0)? 3 : i - 1;

            // Step i RK4 weight.
            w_i = weight[i];

            // Step i time.
            t_i = t_ + w_i;

            // Step i proportion of open synaptic channels.
            Oc_i = Oc + w_i * k_Oc[iprev];

            // Store intermediate results for the KC cell's eventual update.
            gOcs[i] = g * Oc_i;

            // Update RK4 k variable.
            k_Oc[i] = alpha * (1. - Oc_i) * T(t_i) - beta * Oc_i;
        }

        // Update Oc.
        Oc += dt/6 * (k_Oc[0] + 2*k_Oc[1] + 2*k_Oc[2] + k_Oc[3]);
    }

    void CppCholinergicSynapse::generate_random_spikes(double lambda_, unsigned int n_) {
        // Clear any previous values in ts.
        ts.clear();

        // Generate new random activation times if lambda > 0.
        if(lambda > 0) {
            // Handle random number generation.
            std::random_device rd;
            std::mt19937 gen(rd());
            std::exponential_distribution<double> spike_isi(lambda_);

            // Generate random activation times.
            double tprev = 0;
            double tnext;
            for (int i = 0; i < n_; i++) {
                // Compute the next random activation time.
                tnext = tprev + spike_isi(gen);
                // Add t0 first, if it is smaller than the next activation time and larger than the previous.
                if (has_t0 && t0 < tnext && t0 > tprev) {
                    ts.push_back(t0);
                }
                ts.push_back(tnext);
                tprev = tnext;
            }

        // Otherwise, just add in t0.
        } else ts.push_back(t0);

        // Initialize tnow.
        tnow = ts.front();
        ts.pop_front();
    }

    double CppCholinergicSynapse::T(double t_) {
        // Square pulse which models the neurotransmitter opening of the synaptic channels.
        return 0.5 * step(tnow + 0.3 - t_) * step(t_ - tnow);
    }
}
