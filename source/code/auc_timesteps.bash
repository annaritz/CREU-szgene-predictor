set -e  ## this will mean that the script continues


echo FIGURE 4





python3 predict.py   -o ../outfiles/SZ_2-layer_0-sinksource_1500-steps --auc -l 2 --force  --sinksource_method -c 0 -t 2000


python3 predict.py   -o ../outfiles/SZ_2-layer_1-sinksource_1500-steps --auc -l 2 --force  --sinksource_method -c 1 -t 2000






echo DONE.