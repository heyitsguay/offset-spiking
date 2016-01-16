//
// Created by mattguay on 1/8/16.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <random>
#include <cstring>
#include <algorithm>

#include "input_test.h"
#include "time_constant.h"

#define TC 0.63212055882

using namespace kcnet;

double time_constant(double current_A) {
    Network net = input_test(10., 100., current_A, false);
    std::vector<double>::iterator V_max_it = std::max_element(net.Vs.begin(), net.Vs.end());
    // Index of maximum KC membrane potential.
    int V_max_idx = (int)std::distance(net.Vs.begin(), V_max_it);

    // Maximum KC membrane potential.
    double V_max = *V_max_it;

    // Time constant is time from current activation (10ms) until V_tc is
    // reached.
    double V_tc = TC * V_max + (1 - TC) * net.Vs[0];

    // Find where V_tc is reached.
    double d_tc = 99999999.;
    int V_closest_idx = 0;
    for(int i=0; i<V_max_idx; i++) {
        double V = net.Vs[i];
        double d_new = std::abs(V - V_tc);
        if(d_new < d_tc) {
            d_tc = d_new;
            V_closest_idx = i;
        }
    }

    double tc = V_closest_idx * net.dt - 10.;
    std::cout << tc << std::endl;

    return tc;
}

