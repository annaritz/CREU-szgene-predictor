set -e  ## this will mean that the script continues


echo AUC RUNS


python3 predict.py   -o ../outfiles/SZ_2-layer_1-sinksource --auc -l 2 --force -s 50 --sinksource_method -c 1


echo DONE.