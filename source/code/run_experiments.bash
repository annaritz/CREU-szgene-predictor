set -e  ## this will mean that the script continues


echo AUC RUNS


python3 predict.py   -o ../outfiles/SZ_4-layer_0.15-sinksource --auc -l 4 --force -s 50 --sinksource_method -c 0.15 

echo DONE.