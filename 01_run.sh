# Systems
systems=("128w-pos-1")
sizeXs=("15.6404")
sizeYs=("15.6404")
sizeZs=("31.2808")
numAtoms=("384")
# Workflow 
for i in "${!systems[@]}" 
do
    system=${systems[$i]}
    sizeX=${sizeXs[$i]}
    sizeY=${sizeYs[$i]}
    sizeZ=${sizeZs[$i]}
    numAtom=${numAtoms[$i]}

    echo Processing system $system ...
    trajFile=$system.xyz 

#    # Step m1: Resample trajectory 
#    cd m1_resample # (
#    mkdir -p output
#    rm -f $trajFile
#    ln -s ../m2_traj/${trajFile} .
#    for start_frame in 0 10000 20000 30000 40000 50000
#    do
#	subTrajFile=${system}_s${start_frame}.xyz
#	python resample.py $trajFile ${start_frame} output/$subTrajFile
#    done
#    cd .. # ) End Step m1

#    # Step 0: Generate surface trajectory
#    cd 0_prepare # (
#    for start_frame in 0 10000 20000 30000 40000 50000
#    do
#	# a. Prepare input trajectory
#	subTraj=${system}_s${start_frame}
#	subTrajFile=${subTraj}.xyz
#	rm -f $subTrajFile 
#	ln -s ../m1_resample/output/$subTrajFile .
#
#	# b. Prepare input metadata
#	inputFile=input_${subTraj}
#	cp input_template $inputFile
#	sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
#	sed -i "s/SIZEX/${sizeX}/g" $inputFile
#	sed -i "s/SIZEY/${sizeY}/g" $inputFile
#	sed -i "s/SIZEZ/${sizeZ}/g" $inputFile 
#
#	# c. Run Chandler
#	python 1_chandler_fast.py < $inputFile  
#
#	# d. Clean up
#	mkdir -p output
#	mv *.cube output/
#       mv surf_${subTraj}.dat output/  
#    done
#    cd .. # )
#    # End Step 0
    
#    # 1_case1 
#    cd 1_case1 # (
#    # --- Compile and clean up
#    rm -rf ihb_sulpizi 
#    mkdir -p obj
#    make 
#    rm -rf modeule_ihb.mod obj
#    #for start_frame in 0 10000 20000 30000 40000 50000
#    for start_frame in 0 
#    do 
# 		# a. Prepare input trajectory and parameter file
# 		subTraj=${system}_s${start_frame}
# 		subTrajFile=${subTraj}.xyz
# 		surfTrajFile=surf_${subTraj}.dat
# 		rm -rf $subTrajFile $surfTrajFile 
# 		ln -s ../m1_resample/output/$subTrajFile .
# 		ln -s ../0_prepare/output/$surfTrajFile .
# 		for d in {1..6}
# 		do 
# 			inputFile=input_${subTraj}_${d}
# 			cp input_template $inputFile
# 			sed -i "s/SIZEX/${sizeX}/g" $inputFile
# 			sed -i "s/SIZEY/${sizeY}/g" $inputFile
# 			sed -i "s/SIZEZ/${sizeZ}/g" $inputFile
# 			sed -i "s/SYSTEM/${subTraj}/g" $inputFile
# 			sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
# 			sed -i "s/NUMATOM/${numAtom}/g" $inputFile
# 			sed -i "s/SURFTRAJFILE/${surfTrajFile}/g" $inputFile
# 			sed -i "s/THICKNESS/${d}/g" $inputFile
# 		done

# 		# b. Run Chandler
# 		for d in {1..6}
# 		do 
# 			inputFile=input_${subTraj}_${d}
# 			./ihb_sulpizi < $inputFile 
# 		done

# 		# c. Clean up
# 		mkdir -p output
# 		mv ${subTraj}_* output/
# 		mv recentered*.xyz output/
#    done
#    cd .. # ) # End 1_case1

    # 2_case2
    cd 2_case2 # (
    # --- Compile and clean up
    ./compile.sh
    #for start_frame in 0 10000 20000 30000 40000 50000
    for start_frame in 0 
    do 
		# a. Prepare input trajectory and parameter file
		subTraj=${system}_s${start_frame}
		subTrajFile=${subTraj}.xyz
		surfTrajFile=surf_${subTraj}.dat
		rm -rf $subTrajFile $surfTrajFile 
		ln -s ../m1_resample/output/$subTrajFile .
		ln -s ../0_prepare/output/$surfTrajFile .
		for d in {1..6}
		do 
			inputFile=input_${subTraj}_${d}
				cp input_template $inputFile
			sed -i "s/SIZEX/${sizeX}/g" $inputFile
			sed -i "s/SIZEY/${sizeY}/g" $inputFile
			sed -i "s/SIZEZ/${sizeZ}/g" $inputFile
			sed -i "s/SYSTEM/${subTraj}/g" $inputFile
			sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
			sed -i "s/NUMATOM/${numAtom}/g" $inputFile
			sed -i "s/SURFTRAJFILE/${surfTrajFile}/g" $inputFile
			sed -i "s/THICKNESS/${d}/g" $inputFile
		done

		# b. Run ihb 
		for d in {1..6}
		do 
			inputFile=input_${subTraj}_${d}
			./main_interface_c_format2 < $inputFile
			./main_interface_k_format2 < $inputFile
				./main_interface_n_format2 < $inputFile
		done

		# c. Clean up
		mkdir -p output
		mv ${subTraj}_* output/

    done
    cd .. # ) # End 2_case2

    # # 2_orientation
    # cd 2_orientation # (
    # # --- Compile and clean up
    # rm -f iad_sulpizi
	# mkdir -p obj
	# make
	# rm -rf obj
	# rm module_ihb.mod
    # #for start_frame in 0 10000 20000 30000 40000 50000
    # for start_frame in 0 
    # do 
	# 	# a. Prepare input trajectory and parameter file
	# 	subTraj=${system}_s${start_frame}
	# 	subTrajFile=${subTraj}.xyz
	# 	surfTrajFile=surf_${subTraj}.dat
	# 	rm -rf $subTrajFile $surfTrajFile 
	# 	ln -s ../m1_resample/output/$subTrajFile .
	# 	ln -s ../0_prepare/output/$surfTrajFile .
	# 	for d in {1..6}
	# 	do 
	# 		inputFile=input_${subTraj}_${d}
	# 		cp input_template $inputFile
	# 		sed -i "s/SIZEX/${sizeX}/g" $inputFile
	# 		sed -i "s/SIZEY/${sizeY}/g" $inputFile
	# 		sed -i "s/SIZEZ/${sizeZ}/g" $inputFile
	# 		sed -i "s/SYSTEM/${subTraj}/g" $inputFile
	# 		sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
	# 		sed -i "s/NUMATOM/${numAtom}/g" $inputFile
	# 		sed -i "s/SUFTRAJFILE/${surfTrajFile}/g" $inputFile
	# 		sed -i "s/THICKNESS/${d}/g" $inputFile
	# 	done

	# 	# b. Run 
	# 	for d in {1..6}
	# 	do 
	# 		inputFile=input_${subTraj}_${d}
	# 		./iad_sulpizi < $inputFile
	# 	done
		
	# 	# c. Clean up
	# 	mkdir -p output
	# 	mv ${subTraj}_c2*dat output/
	# 	rm gmon.out
	# 	rm recentered_traj_pos_sampled.xyz
	# 	rm ${subTraj}_*_list.dat
    # done
    # cd .. # ) # End 2_case2
done
