cp ../../../build/AMCA3D ./
unzip PBF_x.zip
wait $!
mpirun -np 4 ./AMCA3D -i Inputfile.i -o run.log