set -e  ## this will mean that the script continues


echo FIGURE 4



python3 predict.py   -o ../outfiles/SZ_1-layer_0-sinksource --single -l 1 --force  --sinksource_method -c 0

python3 predict.py   -o ../outfiles/SZ_1-layer_1-sinksource --single -l 1 --force  --sinksource_method -c 1



python3 predict.py   -o ../outfiles/ASD_1-layer_0-sinksource --single -l 1 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/ASD_1-layer_1-sinksource --single -l 1 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt



python3 predict.py   -o ../outfiles/SZ_2-layer_0-sinksource --single -l 2 --force  --sinksource_method -c 0


python3 predict.py   -o ../outfiles/SZ_2-layer_1-sinksource --single -l 2 --force  --sinksource_method -c 1




python3 predict.py   -o ../outfiles/ASD_2-layer_0-sinksource --single -l 2 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt


python3 predict.py   -o ../outfiles/ASD_2-layer_1-sinksource --single -l 2 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt



python3 predict.py   -o ../outfiles/SZ_3-layer_0-sinksource --single -l 3 --force  --sinksource_method -c 0

python3 predict.py   -o ../outfiles/SZ_3-layer_1-sinksource --single -l 3 --force  --sinksource_method -c 1




python3 predict.py   -o ../outfiles/ASD_3-layer_0-sinksource --single -l 3 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt


python3 predict.py   -o ../outfiles/ASD_3-layer_1-sinksource --single -l 3 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

echo DONE.