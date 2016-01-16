//
// Created by mattguay on 1/7/16.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <random>
#include <cstring>
#include <algorithm>

#include "kc_threshold.h"
#include "input_test.h"

using namespace kcnet;

double firing_threshold(Network& net) {

    // Get a sequence of membrane potentials.
    std::vector<double> Vs = net.Vs;
    auto N = (int)Vs.size();

    // Get spike times.
    std::vector<double> spike_times = net.count_spikes();

    if(spike_times.size() > 0) {
        // Index of the first spike in Vs.
        auto idx =  (int)std::round(spike_times[0] / net.dt);

        // Second derivative of Vs (central difference).
        std::vector<double> d2Vs;

        // Calculate second derivative up until the spike.
        for (int i = 1; i < idx; i++) {
            double d2V = Vs[i + 1] - 2 * Vs[i] + Vs[i - 1];
            d2Vs.push_back(d2V);
        }

        // Find the index of the max second derivative before the
        // first spike time.
        std::vector<double>::iterator max = std::max_element(d2Vs.begin(), d2Vs.end());
        int max_idx = (int)std::distance(d2Vs.begin(), max);

        double temp = max_idx * net.dt;

        std::cout << temp << std::endl;

        // Get the membrane potential at that max second derivative index.
        double V_thresh = Vs[max_idx];

        // Get the average membrane potential pre-stimulus.
        double V_avg = 0.;
        // Average values 0 through avg_idx.
        int avg_idx = 300;

        for(int i=0; i<avg_idx; ++i) {
            V_avg += Vs[i];
        }
        V_avg /= (double)avg_idx;

        return V_thresh - V_avg;
    }

}

double kc_threshold(double t0, double t1) {
    Network net = input_test(t0, t1, 0.5);
    double thresh = firing_threshold(net);
    std::cout << thresh << std::endl;
    return thresh;
}




