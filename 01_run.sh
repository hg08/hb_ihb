# Gennerate meta data for trajectory data
python genMeta.py

for system in 128w-pos-1 
do
    echo Processing system $system ...
    trajFile=$system.xyz 
    metaFile=$system.json

    # Step m1: Resample trajectory 
    cd m1_resample # (
    mkdir -p output
    rm -f $trajFile
    ln -s ../m2_traj/${trajFile} .
    for start_frame in 0 10000 20000 30000 40000 50000
    do
	subTrajFile=${system}_s${start_frame}.xyz
	python resample.py $trajFile ${start_frame} output/$subTrajFile
    done
    cd .. # )
    # End Step m1


done

