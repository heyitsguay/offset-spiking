//
// Created by mattguay on 10/28/15.
//
#include "epsp_histogram.h"

//#include <omp.h>
#include <fstream>
#include <sys/time.h>
#include <algorithm>
//#include <cmath>

namespace kcnet {
    void epsp_histogram(const std::string save_name_, const unsigned int n_trials_, const double syn_weight_) {
        /* Generates a PN->KC EPSP distribution histogram with n_trials_ trials, using synaptic weights drawn from an
        empirical distribution weighted by syn_weight_. */

        // Timing variables.
        struct timeval clock0, clock1;
        // Start timing.
        gettimeofday(&clock0, NULL);

        // A single CppNetwork (net) is created with one synapse, and initialized
        // to equilibrium.
        CppNetwork net = CppNetwork();
        double dt = 0.03;
        double t0 = 2;
        double runtime = 15.;
        // Set up net with one synapse.
        net.setup(1, dt, t0, t0, 0., syn_weight_);

        // Load saved equilibrium state.
        net.kc->load_state("eq093015");

        // The peak EPSP is measured as the maximum KC V value following the synapse's activation minus the value just
        // before activation.
        std::vector<double> epsps;

        // Main loop.
        for(int i=0; i<n_trials_; i++) {

            // Run, resetting to the equilibrium state each time.
            net.run(runtime, true);

            // Maximum KC membrane voltage.
            double Vmax = *std::max_element(net.Vs.begin(), net.Vs.end());

            // Membrane voltage right before the synapse activates.
            unsigned int tprev = (unsigned int)std::floor(t0 / dt);
            double Vprev = net.Vs[tprev];

            // Push log10 of EPSP amplitude, expressed in volts.
            epsps.push_back(std::log10((Vmax - Vprev) / 1000.));
        }

        // Write EPSP data to a file.
        std::ofstream write_file;
        std::string fname = "../save/" + save_name_ + ".dat";
        write_file.open(fname.c_str());
        for(int i=0; i<n_trials_; i++) {
            write_file << epsps[i];
            if(i < n_trials_-1) {
                write_file << ",";
            } else write_file << "\n";
        }

        write_file.close();

        // Finish timing information.
        gettimeofday(&clock1, NULL);
        double elapsed = ((clock1.tv_sec - clock0.tv_sec) * 1000000u +
                          clock1.tv_usec - clock0.tv_usec) / 1.e6;

        printf("Wrote %d trials to file %s.dat in %.3f seconds.\n", n_trials_, save_name_.c_str(), elapsed);
    }
}
