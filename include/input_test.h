//
// Created by mattguay on 1/7/16.
//

#ifndef OFFSET_SPIKING_INPUT_TEST_H
#define OFFSET_SPIKING_INPUT_TEST_H
#pragma once

#include "Network.h"

kcnet::Network input_test(
        double current_t0_,
        double current_t1_,
        double current_A_,
        bool plot_flag_=true);

#endif //OFFSET_SPIKING_INPUT_TEST_H
