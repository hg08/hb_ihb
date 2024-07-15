rm -rf ihb_sulpizi *.xyz *.mod *.out *.dat
find . -maxdepth 1 -type f -name 'input_*' ! -name 'input_template' -delete
rm -rf output
