#!/bin/bash

make
./bin/test > data
gnuplot error.gnuplot
