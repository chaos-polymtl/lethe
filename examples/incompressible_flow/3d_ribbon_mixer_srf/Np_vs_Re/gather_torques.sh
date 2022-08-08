#!/bin/bash
for i in $(cat case_index.txt) ; 
do 
echo "Gathering $i"
echo "-------------"
tail -1 ./$i/output/torque.00.dat   >> gather.dat
done