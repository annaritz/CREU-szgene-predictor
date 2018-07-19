set -e  ## this will mean that the script continues


echo FIGURE 2



python3 predict.py   -o ../outfiles/SZ_1-layer_0-sinksource --auc -l 1 --force  --sinksource_method -c 0

python3 predict.py   -o ../outfiles/SZ_1-layer_0.01-sinksource --auc -l 1 --force  --sinksource_method -c 0.01

python3 predict.py   -o ../outfiles/SZ_1-layer_0.1-sinksource --auc -l 1 --force  --sinksource_method -c 0.1

python3 predict.py   -o ../outfiles/SZ_1-layer_1-sinksource --auc -l 1 --force  --sinksource_method -c 1

python3 predict.py   -o ../outfiles/SZ_1-layer_10-sinksource --auc -l 1 --force  --sinksource_method -c 10

python3 predict.py   -o ../outfiles/SZ_1-layer_50-sinksource --auc -l 1 --force  --sinksource_method -c 50



python3 predict.py   -o ../outfiles/ASD_1-layer_0-sinksource --auc -l 1 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_1-layer_0.01-sinksource --auc -l 1 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_1-layer_0.1-sinksource --auc -l 1 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_1-layer_1-sinksource --auc -l 1 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_1-layer_10-sinksource --auc -l 1 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_1-layer_50-sinksource --auc -l 1 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt


echo FIGURE 3

python3 predict.py   -o ../outfiles/SZ_1-layer_0-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0 -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/SZ_1-layer_0.01-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0.01 -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/SZ_1-layer_0.1-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0.1 -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/SZ_1-layer_1-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 1 -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/SZ_1-layer_10-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 10 -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/SZ_1-layer_50-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 50 -w --examine_vs_negatives



python3 predict.py   -o ../outfiles/ASD_1-layer_0-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/ASD_1-layer_0.01-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/ASD_1-layer_0.1-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/ASD_1-layer_1-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/ASD_1-layer_10-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives

python3 predict.py   -o ../outfiles/ASD_1-layer_50-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives


echo FIGURE 6

python3 predict.py   -o ../outfiles/SZ_2-layer_0-sinksource --auc -l 2 --force  --sinksource_method -c 0

python3 predict.py   -o ../outfiles/SZ_2-layer_0.01-sinksource --auc -l 2 --force  --sinksource_method -c 0.01

python3 predict.py   -o ../outfiles/SZ_2-layer_0.1-sinksource --auc -l 2 --force  --sinksource_method -c 0.1

python3 predict.py   -o ../outfiles/SZ_2-layer_1-sinksource --auc -l 2 --force  --sinksource_method -c 1

python3 predict.py   -o ../outfiles/SZ_2-layer_10-sinksource --auc -l 2 --force  --sinksource_method -c 10

python3 predict.py   -o ../outfiles/SZ_2-layer_50-sinksource --auc -l 2 --force  --sinksource_method -c 50



python3 predict.py   -o ../outfiles/ASD_2-layer_0-sinksource --auc -l 2 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_2-layer_0.01-sinksource --auc -l 2 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_2-layer_0.1-sinksource --auc -l 2 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_2-layer_1-sinksource --auc -l 2 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_2-layer_10-sinksource --auc -l 2 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_2-layer_50-sinksource --auc -l 2 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt


python3 predict.py   -o ../outfiles/SZ_3-layer_0-sinksource --auc -l 3 --force  --sinksource_method -c 0

python3 predict.py   -o ../outfiles/SZ_3-layer_0.01-sinksource --auc -l 3 --force  --sinksource_method -c 0.01

python3 predict.py   -o ../outfiles/SZ_3-layer_0.1-sinksource --auc -l 3 --force  --sinksource_method -c 0.1

python3 predict.py   -o ../outfiles/SZ_3-layer_1-sinksource --auc -l 3 --force  --sinksource_method -c 1

python3 predict.py   -o ../outfiles/SZ_3-layer_10-sinksource --auc -l 3 --force  --sinksource_method -c 10

python3 predict.py   -o ../outfiles/SZ_3-layer_50-sinksource --auc -l 3 --force  --sinksource_method -c 50



python3 predict.py   -o ../outfiles/ASD_3-layer_0-sinksource --auc -l 3 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_3-layer_0.01-sinksource --auc -l 3 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_3-layer_0.1-sinksource --auc -l 3 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_3-layer_1-sinksource --auc -l 3 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_3-layer_10-sinksource --auc -l 3 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_3-layer_50-sinksource --auc -l 3 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt


echo DONE.