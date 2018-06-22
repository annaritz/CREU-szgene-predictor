set -e  ## this will mean that the script continues





echo AUCs

python3 predict2.py -d ../infiles/SZ_positives.txt -o ../outfiles/SZ_1-layer_no_neg_0.150 --auc -l 1 --force -s 50 -w

python3 predict2.py -d ../infiles/SZ_positives.txt -o ../outfiles/SZ_1-layer_0.01-sinksource_no_neg_0.150 --auc -l 1 --force -s 50 --sinksource_method -c 0.01 -w

python3 predict2.py -d ../infiles/SZ_positives.txt -o ../outfiles/SZ_1-layer_10-sinksource_no_neg_0.150 --auc -l 1 --force -s 50 --sinksource_method -c 10 -w

python3 predict2.py -d ../infiles/SZ_positives.txt -o ../outfiles/SZ_1-layer_50-sinksource_no_neg_0.150 --auc -l 1 --force -s 50 --sinksource_method -c 50 -w






echo DONE.
