#!/bin/bash
mpicc parallel-nbody.c -std=c11 -lm  -march=native -O3 -o par-nbody