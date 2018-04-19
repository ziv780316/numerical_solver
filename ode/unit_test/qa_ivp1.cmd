#!/bin/bash

debug="-d"
predictor="-p"

# principle test of AB
../bin/ode_solver -m ab -o 1 -t 0.1 ${debug} ${predictor} > output_ivp1_ab1
../bin/ode_solver -m ab -o 2 -t 0.1 ${debug} ${predictor} > output_ivp1_ab2
../bin/ode_solver -m ab -o 3 -t 0.1 ${debug} ${predictor} > output_ivp1_ab3
../bin/ode_solver -m ab -o 4 -t 0.1 ${debug} ${predictor} > output_ivp1_ab4

# stability test of AB
../bin/ode_solver -m ab -o 1 -t 1.9 ${debug} ${predictor} > output_ivp1_ab1_t1p9
../bin/ode_solver -m ab -o 1 -t 2.1 ${debug} ${predictor} > output_ivp1_ab1_t2p1
../bin/ode_solver -m ab -o 2 -t 0.9 ${debug} ${predictor} > output_ivp1_ab2_t0p9
../bin/ode_solver -m ab -o 2 -t 1.1 ${debug} ${predictor} > output_ivp1_ab2_t1p1
../bin/ode_solver -m ab -o 3 -t 0.5 ${debug} ${predictor} > output_ivp1_ab3_t0p5
../bin/ode_solver -m ab -o 3 -t 0.6 ${debug} ${predictor} > output_ivp1_ab3_t0p6

# principle test of AM
../bin/ode_solver -m am -o 1 -t 0.1 ${debug} ${predictor} > output_ivp1_am1
../bin/ode_solver -m am -o 2 -t 0.1 ${debug} ${predictor} > output_ivp1_am2
../bin/ode_solver -m am -o 3 -t 0.1 ${debug} ${predictor} > output_ivp1_am3
../bin/ode_solver -m am -o 4 -t 0.1 ${debug} ${predictor} > output_ivp1_am4

# stability test of AM
../bin/ode_solver -m am -o 1 -t 2.5 ${debug} ${predictor} > output_ivp1_am1_t2p5
../bin/ode_solver -m am -o 2 -t 0.9 ${debug} ${predictor} > output_ivp1_am2_t0p9
../bin/ode_solver -m am -o 2 -t 1.0 ${debug} ${predictor} > output_ivp1_am2_t1p0
../bin/ode_solver -m am -o 2 -t 1.1 ${debug} ${predictor} > output_ivp1_am2_t1p1
../bin/ode_solver -m am -o 3 -t 5.9 ${debug} ${predictor} > output_ivp1_am3_t5p9
../bin/ode_solver -m am -o 3 -t 6.1 ${debug} ${predictor} > output_ivp1_am3_t6p1

# principle test of BDF
../bin/ode_solver -m bdf -o 1 -t 0.1 ${debug} ${predictor} > output_ivp1_bdf1_t0p1
../bin/ode_solver -m bdf -o 2 -t 0.1 ${debug} ${predictor} > output_ivp1_bdf2_t0p1
../bin/ode_solver -m bdf -o 3 -t 0.1 ${debug} ${predictor} > output_ivp1_bdf3_t0p1

# stiff test
../bin/ode_solver -m bdf -o 1 -t 10 ${debug} ${predictor} > output_ivp1_bdf1_t10
../bin/ode_solver -m bdf -o 2 -t 10 ${debug} ${predictor} > output_ivp1_bdf2_t10
../bin/ode_solver -m bdf -o 3 -t 10 ${debug} ${predictor} > output_ivp1_bdf3_t10
../bin/ode_solver -m am -o 2 -t 10 ${debug} ${predictor} > output_ivp1_bdf3_t10
