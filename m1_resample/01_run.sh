original_traj=128w-pos-1.xyz
rm -rf output
mkdir -p output
start_frame=0
sampled_traj=${original_traj%.xyz}-s${start_frame}.xyz
python resample.py ${original_traj} ${start_frame} output/${sampled_traj}

#for start_frame in 0 10000 20000 30000 40000 50000
#do
#    sampled_traj=${original_traj%.xyz}-s${start_frame}.xyz
#    python resample.py ${original_traj} ${start_frame} output/${sampled_traj}
#done

