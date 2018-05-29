#!/bin/bash

#debug="-d"
predictor="-p"

# principle test of RK
../bin/ode_solver -m rk -o 4 -t 0.1 ${debug} ${predictor} > output_ivp4_rk4_t0p1

# stability test of AB
../bin/ode_solver -m ab -o 1 -t 0.001 ${predictor} > output_ivp4_ab1_t0p001
../bin/ode_solver -m ab -o 1 -t 0.1 ${predictor} > output_ivp4_ab1_t0p1
../bin/ode_solver -m ab -o 1 -t 0.5 ${predictor} > output_ivp4_ab1_t0p5
../bin/ode_solver -m ab -o 1 -t 1.0 ${predictor} > output_ivp4_ab1_t1p0
../bin/ode_solver -m ab -o 1 -t 1.9 ${predictor} > output_ivp4_ab1_t1p9
../bin/ode_solver -m ab -o 1 -t 2.0 ${predictor} > output_ivp4_ab1_t2p0
../bin/ode_solver -m ab -o 1 -t 2.1 ${predictor} > output_ivp4_ab1_t2p1

# stability test of AM
../bin/ode_solver -m am -o 1 -t 0.001 ${predictor} > output_ivp4_am1_t0p001
../bin/ode_solver -m am -o 1 -t 0.1 ${predictor} > output_ivp4_am1_t0p1
../bin/ode_solver -m am -o 1 -t 0.5 ${predictor} > output_ivp4_am1_t0p5
../bin/ode_solver -m am -o 1 -t 1.0 ${predictor} > output_ivp4_am1_t1p0
../bin/ode_solver -m am -o 1 -t 1.9 ${predictor} > output_ivp4_am1_t1p9
../bin/ode_solver -m am -o 1 -t 2.0 ${predictor} > output_ivp4_am1_t2p0
../bin/ode_solver -m am -o 0 -t 2.1 ${predictor} > output_ivp4_am1_t2p1

# stability test of BDF
../bin/ode_solver -m bdf -o 2 -t 0.1 ${predictor} > output_ivp4_bdf1_t0p1
../bin/ode_solver -m bdf -o 2 -t 0.5 ${predictor} > output_ivp4_bdf1_t0p5
../bin/ode_solver -m bdf -o 2 -t 1.0 ${predictor} > output_ivp4_bdf1_t1p0
