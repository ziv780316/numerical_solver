#!/bin/bash

debug="-d"
predictor="-p"

# stability test of AB
../bin/ode_solver -m ab -o 1 -t 0.19 ${predictor} > output_ivp3_ab1_t0p19
../bin/ode_solver -m ab -o 1 -t 0.21 ${predictor} > output_ivp3_ab1_t0p21

