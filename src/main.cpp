#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <random>
#include <cstring>

#include "gnuplot-iostream.h"
#include "Network.h"
#include "offset_test.h"
#include "epsp_histogram.h"

using namespace std;
using namespace kcnet;


int main2(double, int, double, double, double, double, double, double);
int main1(std::string, unsigned int, unsigned int, bool);

int main() {
//    main1("memtest", 5, 8, true);
//    main1("test121315_8x_2", 20, 8, true);
//    main1("test121315_8x_3", 20, 8, true);
//    main1("test121315_8x_4", 20, 8, true);
//    main1("test121315_8x_5", 20, 8, true);
//
//    main1("test121515_8x_1_0", 20, 8, false);
//    main1("test121315_8x_2_0", 20, 8, false);
//    main1("test121315_8x_3_0", 20, 8, false);
//    main1("test121315_8x_4_0", 20, 8, false);
//    main1("test121315_8x_5_0", 20, 8, false);

    // lambda=0.0025 --> average spike interarrival time of 400ms --> average of 2.5 spikes/sec at rest.
    main2(100, 200, 0.03, 40., 60., 1.5, 15., 0.0025);
//    main2(100, 200, 0.03, 40., 60., 1.5, 15., 0.);
    return 0;
}

//int main0(int& argc, char& *argv[]) {
//    unsigned int n_trials;
//    double syn_weight;
//    std::string save_name;
//    if(argc == 1) {
//        save_name = "default";
//    } else {
//        save_name = argv[1];
//    }
//    n_trials = (unsigned int) ((argc > 2)? atoi(argv[2]) : 1000);
//    syn_weight = (unsigned int) ((argc > 3)? atof(argv[3]) : 18.);
//    epsp_histogram(save_name, n_trials, syn_weight);
//    return 0;
//}
//
int main1(std::string save_name, unsigned int n_trials, unsigned int n_threads, bool use_noise) {
    /* Runs a series of (n_trials) offset spiking tests using (n_threads) threads.
     * Tests use no background spiking. More test parameters are declared in offset_test.cpp.
     */
    offset_test(save_name, n_trials, n_threads, use_noise);
    return 0;
}

int main2(double runtime, int n_synapses, double dt, double t0, double t1, double sigma_noise, double syn_weight, double lambda) {
    //int n_synapses;
    //double runtime, dt, t0, t1, sigma_noise, syn_weight;
    //bool save_state = true;
//    string state_name = "eq093015";

    random_device rd;
    srand(rd());
    Gnuplot gp;

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
        cout << elapsed << " s elapsed." << endl;

//        if (save_state) { net.kc->save_state(state_name); }

        vector<double> plot_vec0 = net.Vs;
        vector<double> plot_vec1 = net.I_syns;

        vector<pair<double, double> > plot_pts0, plot_pts1;
        for (int i = 0; i < net.ts.size(); i++) {
            double scale;

            plot_pts0.push_back(make_pair(net.ts[i], plot_vec0[i]));

            scale = abs(plot_vec0[i]) > 0.000001? plot_vec0[i] : 0.000001*std::sgn(plot_vec0[i]);
            plot_pts1.push_back(make_pair(net.ts[i], plot_vec1[i]) / scale);
        }

        ostringstream oss0, oss1;

        // Y range boundary points.
        double y0min, y0max, y1min, y1max;
        y0min = *min_element(plot_vec0.begin(), plot_vec0.end()) - 0.0001;
        y0max = *max_element(plot_vec0.begin(), plot_vec0.end()) + 0.0001;
        y1min = *min_element(plot_vec1.begin(), plot_vec1.end()) - 0.0001;
        y1max = *max_element(plot_vec1.begin(), plot_vec1.end()) + 0.0001;

        // Plot command.
        oss0 << "set xrange [0:" << runtime << "]\nset yrange [" << y0min << ":"
        << y0max << "]\n set term qt 0\n";

        oss1 << "set yrange [" << y1min << ":"
        << y1max << "]\n set term qt 1\n";


        gp << oss0.str();

        gp << "plot" << gp.file1d(plot_pts0) << "with lines title 'V'\n";

        gp << oss1.str();

        gp << "plot" << gp.file1d(plot_pts1) << "with lines title 'I_syn'\n";
    }
//    delete(net);

    return 0;
}