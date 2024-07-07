# Systems
systems=("128w-pos-1")
sizeX=("15.6404")
sizeY=("15.6404")
sizeZ=("31.2808")

# Workflow 
for i in "${!systems[@]}" 
do
    system=${systems[$i]}
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

#    # Step 0: Generate surface trajectory
#    cd 0_prepare # (
#    #for start_frame in 0 10000 20000 30000 40000 50000
#    for start_frame in 0 
#    do
#	subTrajFile=${system}_s${start_frame}.xyz
#	rm -f $subTrajFile 
#	ln -s ../m1_resample/output/$subTrajFile .
#	python 1_chandler_fast.py < 
#    done
#    cd .. # )
done

