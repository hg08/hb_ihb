# --- Complile and clear extra files 
rm -rf ihb_sulpizi
mkdir -p obj
make 
rm -rf module_ihb.mod obj

# --- Prepare input files and parameter files
# Trajectory and surface file
rm -rf traj_pos_128w_itp.xyz
rm -rf surf_traj.dat
ln -s ../m1_resample/128w-pos-1_sampled.xyz traj_pos_128w_itp.xyz
ln -s ../0_prepare/surf_128w-pos-1_sampled.dat surf_traj.dat
# Parameter files
for d in {1..6}
do
  cp input_template input_main_sample_128w_itp_$d
  sed -i "s/THICKNESS/$d/g" input_main_sample_128w_itp_$d
  ./ihb_sulpizi < input_main_sample_128w_itp_$d
done

# --- Move output files to output folder
rm -rf output
mkdir -p output
mv 128w_* output/
mv recentered_traj_pos_sampled.xyz output/
