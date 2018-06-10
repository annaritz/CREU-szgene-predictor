set -e  ## this will mean that the script continues

echo SINGLE RUNS

echo SINGLE RUNS
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_Full_Pos_0.100 --single -l 1 --force -s 50
python3 predict2.py -d infiles/SZ_positives_no_iso.txt -o outfiles/SZ_1-layer_no_iso_0.100 --single -l 1 --force -s 50

python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_0.100 --single -l 2 --force -s 50
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_0.100 --single -l 3 --force -s 50

python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_0.1-sinksource_0.100 --single -l 1 --force -s 50 --sinksource_method -c 0.1
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_1-sinksource_0.100 --single -l 1 --force -s 50 --sinksource_method -c 1
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_5-sinksource_0.100 --single -l 1 --force -s 50 --sinksource_method -c 5
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_10-sinksource_0.100 --single -l 1 --force -s 50 --sinksource_method -c 10
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_1-layer_20-sinksource_0.100 --single -l 1 --force -s 50 --sinksource_method -c 20

python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_0.1-sinksource_0.100 --single -l 2 --force -s 50 --sinksource_method -c 0.1
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_1-sinksource_0.100 --single -l 2 --force -s 50 --sinksource_method -c 1
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_5-sinksource_0.100 --single -l 2 --force -s 50 --sinksource_method -c 5
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_10-sinksource_0.100 --single -l 2 --force -s 50 --sinksource_method -c 10
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_20-sinksource_0.100 --single -l 2 --force -s 50 --sinksource_method -c 20

python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_0.1-sinksource_no_neg_0.100 --single -l 2 --force -s 50 --sinksource_method -c 0.1 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_1-sinksource_no_neg_0.100 --single -l 2 --force -s 50 --sinksource_method -c 1 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_5-sinksource_no_neg_0.100 --single -l 2 --force -s 50 --sinksource_method -c 5 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_10-sinksource_no_neg_0.100 --single -l 2 --force -s 50 --sinksource_method -c 10 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_2-layer_20-sinksource_no_neg_0.100 --single -l 2 --force -s 50 --sinksource_method -c 20 -w

python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_0.1-sinksource_0.100 --single -l 3 --force -s 50 --sinksource_method -c 0.1
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_1-sinksource_0.100 --single -l 3 --force -s 50 --sinksource_method -c 1
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_5-sinksource_0.100 --single -l 3 --force -s 50 --sinksource_method -c 5
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_10-sinksource_0.100 --single -l 3 --force -s 50 --sinksource_method -c 10
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_20-sinksource_0.100 --single -l 3 --force -s 50 --sinksource_method -c 20

python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_0.1-sinksource_no_neg_0.100 --single -l 3 --force -s 50 --sinksource_method -c 0.1 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_1-sinksource_no_neg_0.100 --single -l 3 --force -s 50 --sinksource_method -c 1 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_5-sinksource_no_neg_0.100 --single -l 3 --force -s 50 --sinksource_method -c 5 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_10-sinksource_no_neg_0.100 --single -l 3 --force -s 50 --sinksource_method -c 10 -w
python3 predict2.py -d infiles/SZ_positives.txt -o outfiles/SZ_3-layer_20-sinksource_no_neg_0.100 --single -l 3 --force -s 50 --sinksource_method -c 20 -w






echo DONE.