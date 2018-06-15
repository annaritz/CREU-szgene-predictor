set -e  ## this will mean that the script continues

## test on smaller network
python3 predict2.py --single --union --sinksource_constant 10 -l 2 \
	-o combined_outfiles/test_ -g networkfiles/brain_top_geq_0.400.txt \
	-b infiles/motility_positives.txt -d infiles/SZ_positives.txt \
	-n infiles/SZ_negatives.txt --sinksource_method

echo COMBINED RUNS
echo "SZ + CM"
python3 predict2.py --single --union --sinksource_constant 10 -l 2 \
	-o combined_outfiles/SZ_CM -g networkfiles/brain_top_geq_0.150.txt \
	-b infiles/motility_positives.txt -d infiles/SZ_positives.txt \
	-n infiles/SZ_negatives.txt --sinksource_method 

echo "ASD + CM"
python3 predict2.py --single --union --sinksource_constant 10 -l 2 \
	-o combined_outfiles/ASD_CM -g networkfiles/brain_top_geq_0.150.txt \
	-b infiles/motility_positives.txt -d infiles/ASD_positives.txt \
	-n infiles/SZ_negatives.txt --sinksource_method 






