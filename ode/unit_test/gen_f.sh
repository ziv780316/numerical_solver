#!/bin/bash

file="propagate_error_f"
rm -f $file
touch $file
f="exp(-l*x)"
for (( i=1; i<=20; i=i+1 )); do
	f="$f - (d * exp(-l*(x-${i}*h)) * step(x,${i}*h))"
	echo "f${i}(x) = $f" >> $file
done
