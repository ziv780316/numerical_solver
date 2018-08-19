#!/bin/bash
./newton_solver -i normal -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f1_jacobian
./newton_solver -i normal -e forward -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f1_forward
./newton_solver -i normal -e central -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f1_central
./newton_solver -i chord -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f1_chord
./newton_solver -i broyden -e jacobian -m 100 -r 1e-3 -a 1e-6 -u 1e-9 -d -o f1_broyden

for file in `ls f1_*`; do
	diff -q $file golden/$file
done
