#!/bin/bash

debug="-d"
predictor="-p"

# principle test of AB
../bin/ode_solver -m ab -o 1 -t 0.1 ${debug} ${predictor} > output_ivp1_lab1
../bin/ode_solver -m ab -o 2 -t 0.1 ${debug} ${predictor} > output_ivp1_lab2
../bin/ode_solver -m ab -o 3 -t 0.1 ${debug} ${predictor} > output_ivp1_lab3
../bin/ode_solver -m ab -o 4 -t 0.1 ${debug} ${predictor} > output_ivp1_lab4

# stability test of AB
../bin/ode_solver -m ab -o 1 -t 1.9 ${debug} ${predictor} > output_ivp1_lab1_t1p9
../bin/ode_solver -m ab -o 1 -t 2.1 ${debug} ${predictor} > output_ivp1_lab1_t2p1
../bin/ode_solver -m ab -o 2 -t 0.9 ${debug} ${predictor} > output_ivp1_lab2_t0p9
../bin/ode_solver -m ab -o 2 -t 1.1 ${debug} ${predictor} > output_ivp1_lab2_t1p1
../bin/ode_solver -m ab -o 3 -t 0.5 ${debug} ${predictor} > output_ivp1_lab3_t0p5
../bin/ode_solver -m ab -o 3 -t 0.6 ${debug} ${predictor} > output_ivp1_lab3_t0p6

# principle test of AM
../bin/ode_solver -m am -o 1 -t 0.1 ${debug} ${predictor} > output_ivp1_lam1
../bin/ode_solver -m am -o 2 -t 0.1 ${debug} ${predictor} > output_ivp1_lam2
../bin/ode_solver -m am -o 3 -t 0.1 ${debug} ${predictor} > output_ivp1_lam3
../bin/ode_solver -m am -o 4 -t 0.1 ${debug} ${predictor} > output_ivp1_lam4

# stability test of AM
../bin/ode_solver -m am -o 1 -t 2.5 ${debug} ${predictor} > output_ivp1_lam1_t2p5
../bin/ode_solver -m am -o 2 -t 0.9 ${debug} ${predictor} > output_ivp1_lam2_t0p9
../bin/ode_solver -m am -o 2 -t 1.0 ${debug} ${predictor} > output_ivp1_lam2_t1p0
../bin/ode_solver -m am -o 2 -t 1.1 ${debug} ${predictor} > output_ivp1_lam2_t1p1
../bin/ode_solver -m am -o 3 -t 5.9 ${debug} ${predictor} > output_ivp1_lam3_t5p9
../bin/ode_solver -m am -o 3 -t 6.1 ${debug} ${predictor} > output_ivp1_lam3_t6p1
