simulations=("base" "asc" "lb" "asc-lb")

cd data/

for sim in "${simulations[@]}"
do
   echo "Running the $sim simulation"
   time mpirun -np 8 lethe-particles plate-discharge_$sim.prm | tee log2_$sim.out
done
