set -e  ## this will mean that the script continues


echo AUC RUNS


python3 predict.py   -o ../outfiles/SZ_2-layer_0.15-sinksource_hid_negs --roc -l 2 --force  --sinksource_method -c 0.15 

echo DONE.