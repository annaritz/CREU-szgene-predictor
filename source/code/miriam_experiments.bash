set -e  ## this will mean that the script continues



echo FIGURE 6 ASD



python3 predict.py   -o ../outfiles/new_runs/ASD_2-layer_0-sinksource --auc -l 2 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/new_runs/ASD_2-layer_0.01-sinksource --auc -l 2 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/new_runs/ASD_2-layer_0.1-sinksource --auc -l 2 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/new_runs/ASD_2-layer_1-sinksource --auc -l 2 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/new_runs/ASD_2-layer_10-sinksource --auc -l 2 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

python3 predict.py   -o ../outfiles/new_runs/ASD_2-layer_50-sinksource --auc -l 2 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt



echo DONE.