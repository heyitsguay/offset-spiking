#include <iostream>

#include "clamp_test.h"
#include "epsp_histogram.h"
#include "run_once.h"
#include "offset_test.h"
#include "input_test.h"
#include "kc_threshold.h"
#include "time_constant.h"

using namespace kcnet;

int main() {
//    offset_test("memtest", 5, 8, true);
//    offset_test("test121315_8x_2", 20, 8, true);
//    offset_test("test121315_8x_3", 20, 8, true);
//    offset_test("test121315_8x_4", 20, 8, true);
//    offset_test("test121315_8x_5", 20, 8, true);
//
//    offset_test("test121515_8x_1_0", 20, 8, false);
//    offset_test("test121315_8x_2_0", 20, 8, false);
//    offset_test("test121315_8x_3_0", 20, 8, false);
//    offset_test("test121315_8x_4_0", 20, 8, false);
//    offset_test("test121315_8x_5_0", 20, 8, false);

//    run_once(100, 200, 0.03, 40., 60., 1.5, 15., 0.0025);

//    input_test(10., 50., 0.4);

//    kc_threshold();

//    time_constant(0.025);

    clamp_test(-10., 100., 0.003, 1., 0.3, 0.1);

    return 0;
}