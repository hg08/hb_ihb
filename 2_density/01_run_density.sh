rm -rf 128w-pos-1.xyz density
ln -s ../m2_traj/128w-pos-1.xyz .

#Compile
gfortran -o density density.f95


# Run
./density < input_density_OH 
