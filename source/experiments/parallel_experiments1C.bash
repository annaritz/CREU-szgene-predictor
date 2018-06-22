set -e  ## this will mean that the script continues





echo AUCs


python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_0.01-sinksource_0.150 --auc -l 3 --force -s 50 --sinksource_method -c 0.01
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_0.01-sinksource_0.150 --auc -l 1 --force -s 50 --sinksource_method -c 0.01 -w




echo DONE.
