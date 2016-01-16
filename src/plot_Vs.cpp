//
// Created by mattguay on 1/10/16.
//

#include <cstring>
#include <vector>
#include <sstream>


#include "gnuplot-iostream.h"
#include "plot_Vs.h"

// Signum function.
template <typename T> int sgn(T x) {
    return (T(0) < x) - (x < T(0));
}

namespace kcnet {

    void plot_Vs(Network& net) {
        Gnuplot gp;

        std::vector<double> plot_vec0 = net.Vs;

        // Runtime that produced those Vs.
        double runtime = net.Vs.size() * net.dt;

        std::vector<std::pair<double, double> > plot_pts0;
        for (int i = 0; i < net.ts.size(); i++) {
            plot_pts0.push_back(std::make_pair(net.ts[i], plot_vec0[i]));
        }

        std::ostringstream oss0;

        // Y range boundary points.
        double y0min, y0max;
        y0min = *min_element(plot_vec0.begin(), plot_vec0.end()) - 0.0001;
        y0max = *max_element(plot_vec0.begin(), plot_vec0.end()) + 0.0001;

        // Plot command.
        oss0 << "set xrange [0:" << runtime << "]\nset yrange [" << y0min << ":"
        << y0max << "]\n set term qt 0\n";

        gp << oss0.str();

        gp << "plot" << gp.file1d(plot_pts0) << "with lines title 'V'\n";
    }

}