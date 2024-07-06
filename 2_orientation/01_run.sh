# 1. Compile
rm -f iad_sulpizi
mkdir -p obj
make
rm -rf obj
rm module_ihb.mod



# 2. Prepare input
rm -f traj_pos.xyz surf_traj.dat
ln -s ../m1_resample/128w-pos-1_sampled.xyz traj_pos.xyz
ln -s ../0_prepare/surf_128w-pos-1_sampled.dat surf_traj.dat

for d in {1..6}
do
  cp input_relax_X  input_relax_${d}
  sed -i "s/THICKNESS/$d/g" input_relax_${d}
done

# 3. Run
for d in {1..6}
do
  ./iad_sulpizi < input_relax_${d}
done

# 4. Clean
mkdir -p output
mv 128w_itp_c2*dat output
rm input_relax_[1-6]
rm gmon.out
rm recentered_traj_pos_sampled.xyz 
rm 128w_itp_OH_at_interface_list.dat


