#!/bin/bash
for i in $(cat case_index.txt) ; 
do 
echo "Gathering $i"
echo "-------------"
tail -1 ./$i/output/torque.01.dat   >> gather.dat
done