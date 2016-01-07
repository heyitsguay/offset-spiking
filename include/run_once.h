//
// Created by mattguay on 1/7/16.
//

#ifndef OFFSET_SPIKING_RUN_ONCE_H
#define OFFSET_SPIKING_RUN_ONCE_H
#pragma once

int run_once(double runtime,
             int n_synapses,
             double dt,
             double t0,
             double t1,
             double sigma_noise,
             double syn_weight,
             double lambda);


#endif //OFFSET_SPIKING_RUN_ONCE_H
