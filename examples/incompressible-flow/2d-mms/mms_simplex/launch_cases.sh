 for i in $(ls . | grep mms_2d_steady_mesh_ref); do echo $i ; cd $i ; mpirun -np $1 lethe-fluid mms_2d_steady.prm ; cd ..;  done 
