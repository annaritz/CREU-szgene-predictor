set -e  ## this will mean that the script continues

## test on smaller network
python3 predict2.py --single --union --sinksource_constant 10 -l 2 \
	-o ../outfiles/combined_test_ -g ../infiles/networkfiles/brain_top_geq_0.400.txt \
	-b ../infiles/motility_positives.txt -d ../infiles/SZ_positives.txt \
	-n ../infiles/SZ_negatives.txt --sinksource_method

echo COMBINED RUNS
echo "SZ + CM"
python3 predict2.py --single --union --sinksource_constant 10 -l 2 \
	-o ../outfiles/combined_SZ_CM -g ../infiles/networkfiles/brain_top_geq_0.150.txt \
	-b ../infiles/motility_positives.txt -d ../infiles/SZ_positives.txt \
	-n ../infiles/SZ_negatives.txt --sinksource_method 

echo "ASD + CM"
python3 predict2.py --single --union --sinksource_constant 10 -l 2 \
	-o ../outfiles/combined_ASD_CM -g ../infiles/networkfiles/brain_top_geq_0.150.txt \
	-b ../infiles/motility_positives.txt -d ../infiles/ASD_positives.txt \
	-n ../infiles/SZ_negatives.txt --sinksource_method 






