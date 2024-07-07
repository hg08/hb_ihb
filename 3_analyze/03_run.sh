mkdir -p output
rm -f *png

echo "Processing C2 decay rate..."

# Fit decay constant tau2 of C2 for all layers
tau2=output/tau2.dat
rm -rf $tau2
for i in {1..6}
do
	filePath=../2_orientation/output
	./fit_c2.py ${filePath}/128w_itp_c2_ihb_${i}*.dat >> $tau2
done

# Add the THICKNESS
awk '{print NR, $0}' $tau2 > tmp
mv tmp $tau2

