#!/bin/bash

echo `hostname`


list="A14
A7
C0
C14
C7
B14
B7"

for file in $list; do
        sbatch $1 $file
done

echo "=========== $1 Job Submitted  ================="


