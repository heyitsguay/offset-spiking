//
// Created by Matthew Guay on 10/11/15.
//
#pragma once
#include "network_cpp.h"

namespace kcnet {
    void offset_test(std::string save_name_, const unsigned int n_trials_=10, const unsigned int n_threads_=1, const bool use_noise_=true);
}