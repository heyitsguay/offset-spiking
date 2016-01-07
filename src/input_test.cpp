//
// Created by mattguay on 1/7/16.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <random>
#include <cstring>

#include "input_test.h"
#include "gnuplot-iostream.h"
#include "Network.h"

using namespace kcnet;

void input_test(double current_t0_, double current_t1_, double current_A_) {
    Gnuplot gp;

    // Create a Network.
    auto net = Network();

    // Network parameters.
    int n_synapses = 0;

    // RK4 time step.
    double dt = 0.03;

    // Synapse activation window endpoints - not doing anything in this case.
    double t0 = 0.;
    double t1 = 1.;

    // Network setup.
    net.setup(n_synapses, dt, t0, t1);

    // Add input current.
    net.input_current(current_t0_, current_t1_, current_A_);

    // Load KC equilibrium state.
    net.kc->load_state("equilibrium");

    // Network runtime.
    double runtime = current_t1_ + 100.;

    // Run.
    net.run(runtime, true);

    // Plot KC membrane potential.
    std::vector<double> plot_vec0 = net.Vs;

    // (x,y) pairs of points to plot.
    std::vector<std::pair<double, double>> plot_pts;

    // Create (x,y) pairs.
    for (int i = 0; i < net.ts.size(); i++) {
        plot_pts.push_back(std::make_pair(net.ts[i], plot_vec0[i]));
    }

    // String stream to push to GNUPlot.
    std::ostringstream oss;

    // Y range boundary points.
    double y0min, y0max;
    y0min = *min_element(plot_vec0.begin(), plot_vec0.end()) - 0.0001;
    y0max = *max_element(plot_vec0.begin(), plot_vec0.end()) + 0.0001;

    // Plot command string.
    oss << "set xrange [0:" << runtime << "]\nset yrange [" << y0min << ":"
    << y0max << "]\n set term qt 0\n";

    // Send command to GNUplot.
    gp << oss.str();

    gp << "plot" << gp.file1d(plot_pts) << "with lines title 'V'\n";
}