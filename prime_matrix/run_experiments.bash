set -e  ## this will mean that the script continues

echo SINGLE RUNS


echo AUCs

python3 predict2.py -o outfiles/2-layer_20_sinksource_with_negatives --auc --sinksource_constant 10 -l 2 --force --sinksource_method -s 15


echo DONE.