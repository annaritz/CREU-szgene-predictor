set -e  ## this will mean that the script continues





echo singles




python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_0.01-sinksource_0.150 --single -l 2 --force  --sinksource_method -c 0.01

python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_0.01-sinksource_0.150 --single -l 1 --force  --sinksource_method -c 0.01 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_10-sinksource_0.150 --single -l 1 --force  --sinksource_method -c 10 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_50-sinksource_0.150 --single -l 1 --force  --sinksource_method -c 50 -w

python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_0.01-sinksource_no_neg_0.150 --single -l 2 --force  --sinksource_method -c 0.01 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_10-sinksource_no_neg_0.150 --single -l 2 --force  --sinksource_method -c 10 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_50-sinksource_no_neg_0.150 --single -l 2 --force  --sinksource_method -c 50 -w

python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_0.01-sinksource_no_neg_0.150 --single -l 3 --force  --sinksource_method -c 0.01 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_10-sinksource_no_neg_0.150 --single -l 3 --force  --sinksource_method -c 10 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_50-sinksource_no_neg_0.150 --single -l 3 --force  --sinksource_method -c 50 -w




echo DONE.
