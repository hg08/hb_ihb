# --- Compile four programs rm -rf main_interface_nf_format2
rm -rf main_interface_nf_format2
cd src
# 1. nf
gfortran -Wall -fcheck=all module_tools.f95 types.f95 \
  parameter_shared.f95 atom_module.f95 surf_module.f95 \
  water_module.f95 water_pair_module.f95 count_time.f95 \
  traj_format2.f95 surf_traj.f95 read_interface_input.f95 \
  load.f95 ghbond.f95 gfreeoh.f95 \
  ghbacf_interface_nf_pbc_format2.f95 main_interface_nf_format2.f95 \
  -o main_interface_nf_format2


# clean and move
rm *.mod
mv main_interface_nf_format2 ../
cd ..
