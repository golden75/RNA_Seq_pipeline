#!/bin/bash

echo `hostname`


list="Uribe_014
Uribe_07
Uribe_C0
Uribe_C14
Uribe_C7
Uribe_P14
Uribe_P7"

for file in $list; do
        sbatch $1 $file
done

echo "=========== $1 Job Submitted  ================="


