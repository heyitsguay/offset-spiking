//
// Created by mattguay on 1/10/16.
//

#include "Network.h"
#include "clamp_test.h"
#include "plot_Vs.h"

using namespace kcnet;

void clamp_test(
        double V_com_,
        double runtime_,
        double dt_,
        double sigma_noise_,
        double t_i_,
        double t_d_
        ) {

    // Create a Network.
    auto net = Network();
    // Setup.
    net.setup(0, dt_, 0., 0., sigma_noise_);
    // Set the KC to equilibrium.
    net.kc->load_state("equilibrium");
    // Add a VoltageClamp.
    net.add_VoltageClamp(V_com_, t_i_, t_d_);
    // Set VoltageClamp PID type.
    net.vclamp->pid_type = PD_PV;

    // Run the Network.
    net.run(runtime_, true);

    // Plot the KC membrane potential.
    plot_Vs(net);
}
