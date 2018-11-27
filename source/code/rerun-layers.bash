set -e  ## this will mean that the script continues

echo =========== AUC 50 ==========
python3 -u predict.py -o ../outfiles/ASD_3-layer_50-sinksource --single --union -l 3 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt


echo =========== AUC 10 ==========
python3 -u predict.py -o ../outfiles/ASD_3-layer_10-sinksource --single --union -l 3 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt


echo =========== AUC 1 ==========
python3 -u predict.py -o ../outfiles/ASD_3-layer_1-sinksource --single --union -l 3 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt