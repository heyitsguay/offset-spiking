//
// Created by mattguay on 1/7/16.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <random>
#include <cstring>

#include "plot_Vs.h"
#include "Network.h"
#include "run_once.h"

using namespace kcnet;

void run_once(double runtime, int n_synapses, double dt, double t0, double t1, double sigma_noise, double syn_weight, double lambda) {

    Network net = Network();
    net.setup(n_synapses, dt, t0, t1, sigma_noise, syn_weight, lambda);
    net.create_synapses(400 - n_synapses, syn_weight, false);

    net.kc->load_state("equilibrium");

    clock_t ct0, ct1;
    double elapsed;
    for(int j=0; j<1; j++) {
        ct0 = clock();
        net.run(runtime, true);
        ct1 = clock();
        elapsed = (ct1 - ct0) / 1000000.;
        std::cout << elapsed << " s elapsed." << std::endl;
    }

    // Plot KC membrane potential.
    plot_Vs(net);
}
