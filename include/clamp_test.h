//
// Created by mattguay on 1/10/16.
//

#ifndef OFFSET_SPIKING_CLAMP_TEST_H
#define OFFSET_SPIKING_CLAMP_TEST_H
#pragma once

void clamp_test(
        double V_com_,
        double runtime_,
        double dt_=0.03,
        double sigma_noise_=1.,
        double t_i_=-1.,
        double t_d_=0.
        );

#endif //OFFSET_SPIKING_CLAMP_TEST_H
