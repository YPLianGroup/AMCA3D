cp ../../../build/AMCA3D ./
unzip PBF_y.zip
wait $!
mkdir PreRun
cd PreRun
cp ../AMCA3D ./
mpirun -np 4 ./AMCA3D -i ../Inputfile_0.i -o run_0.log
wait $!
cd ..
mpirun -np 4 ./AMCA3D -i Inputfile_1.i -o run_1.log