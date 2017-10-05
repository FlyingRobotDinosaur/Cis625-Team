#!/bin/bash
let i=10
while [ $i -lt 999999 ]; do
mpirun -np 8 /homes/aan/hw3/MPI_keys 50000 $i 8 8
let i=i*10
done
