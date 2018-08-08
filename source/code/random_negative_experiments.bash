set -e  ## this will mean that the script continues
## the python3 -u "unbuffers" the python call so it flushes to output file.
## note that there are NO --force arguments here.
# Usage:
# bash random_negative_experiments.bash > random_negative_experiments.out

echo RANDOM NEGATIVES

echo SZ

python3 -u predict.py   -o ../outfiles/SZ_1-layer_0-sinksource-random_negatives --auc -l 1  --sinksource_method -c 0 --random_negatives
python3 -u predict.py   -o ../outfiles/SZ_1-layer_0.01-sinksource-random_negatives --auc -l 1  --sinksource_method -c 0.01 --random_negatives
python3 -u predict.py   -o ../outfiles/SZ_1-layer_0.1-sinksource-random_negatives --auc -l 1  --sinksource_method -c 0.1 --random_negatives
python3 -u predict.py   -o ../outfiles/SZ_1-layer_1-sinksource-random_negatives --auc -l 1  --sinksource_method -c 1 --random_negatives
python3 -u predict.py   -o ../outfiles/SZ_1-layer_10-sinksource-random_negatives --auc -l 1  --sinksource_method -c 10 --random_negatives
python3 -u predict.py   -o ../outfiles/SZ_1-layer_50-sinksource-random_negatives --auc -l 1  --sinksource_method -c 50 --random_negatives

echo ASD

python3 -u predict.py   -o ../outfiles/ASD_1-layer_0-sinksource-random_negatives --auc -l 1  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py   -o ../outfiles/ASD_1-layer_0.01-sinksource-random_negatives --auc -l 1  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py   -o ../outfiles/ASD_1-layer_0.1-sinksource-random_negatives --auc -l 1  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py   -o ../outfiles/ASD_1-layer_1-sinksource-random_negatives --auc -l 1  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py   -o ../outfiles/ASD_1-layer_10-sinksource-random_negatives --auc -l 1  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py   -o ../outfiles/ASD_1-layer_50-sinksource-random_negatives --auc -l 1  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives

echo DEGREE PRESERVING RANDOM NEGATIVES

echo SZ

python3 -u predict.py   -o ../outfiles/SZ_1-layer_0-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 0 --random_negatives_degree
python3 -u predict.py   -o ../outfiles/SZ_1-layer_0.01-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 0.01 --random_negatives_degree
python3 -u predict.py   -o ../outfiles/SZ_1-layer_0.1-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 0.1 --random_negatives_degree
python3 -u predict.py   -o ../outfiles/SZ_1-layer_1-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 1 --random_negatives_degree
python3 -u predict.py   -o ../outfiles/SZ_1-layer_10-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 10 --random_negatives_degree
python3 -u predict.py   -o ../outfiles/SZ_1-layer_50-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 50 --random_negatives_degree

echo ASD

python3 -u predict.py   -o ../outfiles/ASD_1-layer_0-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py   -o ../outfiles/ASD_1-layer_0.01-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py   -o ../outfiles/ASD_1-layer_0.1-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py   -o ../outfiles/ASD_1-layer_1-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py   -o ../outfiles/ASD_1-layer_10-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py   -o ../outfiles/ASD_1-layer_50-sinksource-random_negatives_degree --auc -l 1  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree


echo DONE.