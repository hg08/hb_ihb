for folder in m1_resample 0_prepare 1_case1 2_case2 2_orientation 3_analyze 4_plot
do
    cd $folder
    ./clean.sh
    cd ..
done

rm -rf *.out *.err
