//
// Created by mattguay on 1/8/16.
//

#ifndef OFFSET_SPIKING_VOLTAGECLAMP_H
#define OFFSET_SPIKING_VOLTAGECLAMP_H
#pragma once

#include "KenyonCell.h"

namespace kcnet {

    // Dictates the possible PID types.
    enum PID {
        NO_PV, // No process variable (PV) terms.
        D_PV, // PV term on derivative.
        PD_PV // PV term on derivative and proportion.
    };

    class VoltageClamp {
    public:
        PID pid_type; // Current PID type

        double V_final; // Final command potential.
        double V_com; // Command potential at the current time step.
        double I_clamp; // Clamp current.
        double t_i; // Integral time.
        double t_d; // Derivative time.
        double K_p; // Proportional gain.
        double e; // error[n] - current step error.
        double tau; // low pass filter weight.

        VoltageClamp(double V_final_, KenyonCell& kc_, double t_i_=-1, double t_d_=0.);

        ~VoltageClamp();

        void update();

    protected:
        KenyonCell& kc; // Network this thing is interacting with.
        double& V_in; // KC membrane potential.
        double V_oldfinal; // temporary, hopefully...
        double& dt; // RK4 time step.
        double e1; // error[n-1]
        double e2; // error[n-2]
        double td_dt; // Internal variable for performance.
        double dt_ti; // Internal variable for performance.
        double tau_dt; // Internal variable for performance.
        double v1; // Internal variable used when PID has PVs.
        double expt; // RC filter weight.
    };
}

#endif //OFFSET_SPIKING_VOLTAGECLAMP_H
