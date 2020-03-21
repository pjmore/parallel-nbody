#!/bin/bash
if [ -z "$1" ]
then
P=4
else
P=$1
fi


mpirun -np $P ./par-nbody './sample-input.txt' 0.1 199584000 500000