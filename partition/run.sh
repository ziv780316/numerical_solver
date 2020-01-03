#/bin/bash

all_problem=`ls *.so`
for problem in $all_problem; do
	./newton_solver -p ./$problem -e forward -i normal -d -o ${problem/.so/.log}
done
