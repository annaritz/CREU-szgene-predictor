set -e  ## this will mean that the script continues


echo AUC RUNS

python3 predict.py   -o ../outfiles/SZ_4-layer_0.15-sinksource_0.400E --auc -l 4 --force -s 2 --sinksource_method -c 0.15 -d ../infiles/SZ_positives_E.txt -b ../infiles/motility_positives_E.txt -n ../infiles/SZ_negatives_E.txt







echo DONE.