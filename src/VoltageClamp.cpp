//
// Created by mattguay on 1/8/16.
//
#include <cmath>
#include <cassert>

#include "VoltageClamp.h"

namespace kcnet {

    VoltageClamp::VoltageClamp(
            double V_final_,
            KenyonCell& kc_,
            double t_i_,
            double t_d_
            ) :
        V_final(V_final_), kc(kc_), V_in(kc.V), dt(kc.dt) {

        assert(t_d_ >= 0);

        V_oldfinal = kc.V;

        // Initialize clamp current.
        I_clamp = 0.;

        // Initialize command potential and error terms.
        V_com = e = e1 = e2 = 0.;

        // Set PID type.
        pid_type = NO_PV;

        // Initialize clamp control constants.
        t_i = (t_i_ < 0)? dt : t_i_; // Default is dt.
        t_d = t_d_;
        K_p = kc.Cm / dt;
        tau = 5 * dt;
        expt = (dt / tau > 1e-15)? std::exp(-dt / tau) : 1 - dt / tau;

        // Set these so that we don't have to keep doing divisions.
        dt_ti = dt / t_i;
        td_dt = t_d / dt;
        tau_dt = tau / dt;
    }

    VoltageClamp::~VoltageClamp() {}

    void VoltageClamp::update() {
    // Update the VoltageClamp states one time step.

        double d_cmd = V_final - V_oldfinal;
        V_com = V_final + d_cmd * (1 - tau_dt) + (V_com - V_final + d_cmd * tau_dt) * expt;
        V_oldfinal = V_final;

        // Difference between command voltage and membrane potential.
        e = V_com - V_in;

        if(pid_type == NO_PV) {
            I_clamp += K_p * ((1. + dt_ti + td_dt) * e
                       - (1. + 2. * td_dt) * e1 + td_dt * e2);
            e2 = e1;
            e1 = e;

        } else if(pid_type == D_PV) {
            I_clamp += K_p * ((1. + dt_ti) * e - e1
                       + td_dt * (V_in - 2 * v1 + e2));
            e2 = v1;
            v1 = V_in;
            e1 = e;

        } else if(pid_type == PD_PV) {
            I_clamp += K_p * (V_in - v1 + dt_ti * e
                       + td_dt * (V_in - 2 * v1 + e2));
            e2 = v1;
            v1 = V_in;
        }
    }
}
