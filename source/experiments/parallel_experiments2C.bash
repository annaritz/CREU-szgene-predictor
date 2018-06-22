set -e  ## this will mean that the script continues





echo AUCs

python3 predict2.py -d ../infiles/SZ_positives.txt -o ../outfiles/SZ_3-layer_10-sinksource_0.150 --auc -l 3 --force -s 50 --sinksource_method -c 10
python3 predict2.py -d ../infiles/SZ_positives.txt -o ../outfiles/SZ_1-layer_10-sinksource_0.150 --auc -l 1 --force -s 50 --sinksource_method -c 10 -w







echo DONE.
