set -e  ## this will mean that the script continues

echo SINGLE RUNS

echo AUCs
python3 predict2.py -o outfiles/2-layer_5_sinksource --auc --sinksource_constant 5 -l 2 --force
python3 predict2.py -o outfiles/2-layer_10_sinksource --auc --sinksource_constant 10 -l 2 --force
python3 predict2.py -o outfiles/2-layer_20_sinksource --auc --sinksource_constant 20 -l 2 --force
python3 predict2.py -o outfiles/3-layer_10_sinksource --auc --sinksource_constant 10 -l 3 --force




## aggregate
# python3 predict2.py --aggregate outfiles/2-layer_5_sinksource --aggregate outfiles/2-layer_10_sinksource --aggregate outfiles/2-layer_20_sinksource --aggregate outfiles/3-layer_10_sinksource --aggregate_names 2Five --aggregate_names 2Ten --aggregate_names 2Twenty --aggregate_names 3Ten
echo DONE.
