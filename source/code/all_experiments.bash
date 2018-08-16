set -e  ## this will mean that the script continues
# the -u spits output to file without buffering (so you can see where each program is)


## LAYER 1 AUC
python3 -u predict.py -o ../outfiles/SZ_1-layer_0-sinksource --auc -l 1 --force  --sinksource_method -c 0
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.01-sinksource --auc -l 1 --force  --sinksource_method -c 0.01
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.1-sinksource --auc -l 1 --force  --sinksource_method -c 0.1
python3 -u predict.py -o ../outfiles/SZ_1-layer_1-sinksource --auc -l 1 --force  --sinksource_method -c 1
python3 -u predict.py -o ../outfiles/SZ_1-layer_10-sinksource --auc -l 1 --force  --sinksource_method -c 10
python3 -u predict.py -o ../outfiles/SZ_1-layer_50-sinksource --auc -l 1 --force  --sinksource_method -c 50

python3 -u predict.py -o ../outfiles/ASD_1-layer_0-sinksource --auc -l 1 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.01-sinksource --auc -l 1 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.1-sinksource --auc -l 1 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_1-sinksource --auc -l 1 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_10-sinksource --auc -l 1 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_50-sinksource --auc -l 1 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

## LAYER 1 SS+ AUC (no negatives)
python3 -u predict.py -o ../outfiles/SZ_1-layer_0-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0 -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.01-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0.01 -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.1-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0.1 -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_1-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 1 -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_10-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 10 -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_50-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 50 -w --examine_vs_negatives

python3 -u predict.py -o ../outfiles/ASD_1-layer_0-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.01-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.1-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_1-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_10-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_50-sinksource_no_neg --auc -l 1 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt -w --examine_vs_negatives

## LAYER 1 RAND NEGS AUC
python3 -u predict.py -o ../outfiles/SZ_1-layer_0-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 0 --random_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.01-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 0.01 --random_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.1-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 0.1 --random_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_1-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 1 --random_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_10-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 10 --random_negatives
python3 -u predict.py -o ../outfiles/SZ_1-layer_50-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 50 --random_negatives

python3 -u predict.py -o ../outfiles/ASD_1-layer_0-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.01-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.1-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_1-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_10-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives
python3 -u predict.py -o ../outfiles/ASD_1-layer_50-sinksource-random_negatives --auc -l 1 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives

## LAYER 1 RAND NEGS PRESERVING DEGREE AUC
python3 -u predict.py -o ../outfiles/SZ_1-layer_0-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 0 --random_negatives_degree
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.01-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 0.01 --random_negatives_degree
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.1-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 0.1 --random_negatives_degree
python3 -u predict.py -o ../outfiles/SZ_1-layer_1-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 1 --random_negatives_degree
python3 -u predict.py -o ../outfiles/SZ_1-layer_10-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 10 --random_negatives_degree
python3 -u predict.py -o ../outfiles/SZ_1-layer_50-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 50 --random_negatives_degree

python3 -u predict.py -o ../outfiles/ASD_1-layer_0-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.01-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.1-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py -o ../outfiles/ASD_1-layer_1-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py -o ../outfiles/ASD_1-layer_10-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree
python3 -u predict.py -o ../outfiles/ASD_1-layer_50-sinksource-random_negatives_degree --auc -l 1 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt --random_negatives_degree

## LAYER 1 OLD NEGS AUC
python3 -u predict.py -o ../outfiles/SZ_1-layer_0-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 0 -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.01-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 0.01 -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.1-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 0.1 -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/SZ_1-layer_1-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 1 -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/SZ_1-layer_10-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 10 -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/SZ_1-layer_50-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 50 -n ../infiles/old_negatives.txt

python3 -u predict.py -o ../outfiles/ASD_1-layer_0-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.01-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.1-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_1-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_10-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/old_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_50-sinksource-old_negatives --auc -l 1 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/old_negatives.txt

## LAYER 2 AUC
python3 -u predict.py -o ../outfiles/SZ_2-layer_0-sinksource --auc -l 2 --force  --sinksource_method -c 0
python3 -u predict.py -o ../outfiles/SZ_2-layer_0.01-sinksource --auc -l 2 --force  --sinksource_method -c 0.01
python3 -u predict.py -o ../outfiles/SZ_2-layer_0.1-sinksource --auc -l 2 --force  --sinksource_method -c 0.1
python3 -u predict.py -o ../outfiles/SZ_2-layer_1-sinksource --auc -l 2 --force  --sinksource_method -c 1
python3 -u predict.py -o ../outfiles/SZ_2-layer_10-sinksource --auc -l 2 --force  --sinksource_method -c 10
python3 -u predict.py -o ../outfiles/SZ_2-layer_50-sinksource --auc -l 2 --force  --sinksource_method -c 50

python3 -u predict.py -o ../outfiles/ASD_2-layer_0-sinksource --auc -l 2 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_0.01-sinksource --auc -l 2 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_0.1-sinksource --auc -l 2 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_1-sinksource --auc -l 2 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_10-sinksource --auc -l 2 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_50-sinksource --auc -l 2 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

## LAYER 3 AUC
python3 -u predict.py -o ../outfiles/SZ_3-layer_0-sinksource --auc -l 3 --force  --sinksource_method -c 0
python3 -u predict.py -o ../outfiles/SZ_3-layer_0.01-sinksource --auc -l 3 --force  --sinksource_method -c 0.01
python3 -u predict.py -o ../outfiles/SZ_3-layer_0.1-sinksource --auc -l 3 --force  --sinksource_method -c 0.1
python3 -u predict.py -o ../outfiles/SZ_3-layer_1-sinksource --auc -l 3 --force  --sinksource_method -c 1
python3 -u predict.py -o ../outfiles/SZ_3-layer_10-sinksource --auc -l 3 --force  --sinksource_method -c 10
python3 -u predict.py -o ../outfiles/SZ_3-layer_50-sinksource --auc -l 3 --force  --sinksource_method -c 50

python3 -u predict.py -o ../outfiles/ASD_3-layer_0-sinksource --auc -l 3 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_0.01-sinksource --auc -l 3 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_0.1-sinksource --auc -l 3 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_1-sinksource --auc -l 3 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_10-sinksource --auc -l 3 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_50-sinksource --auc -l 3 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

## LAYER 1 SINGLE RUNS
python3 -u predict.py -o ../outfiles/SZ_1-layer_0-sinksource --single --union -l 1 --force  --sinksource_method -c 0
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.01-sinksource --single --union -l 1 --force  --sinksource_method -c 0.01
python3 -u predict.py -o ../outfiles/SZ_1-layer_0.1-sinksource --single --union -l 1 --force  --sinksource_method -c 0.1
python3 -u predict.py -o ../outfiles/SZ_1-layer_1-sinksource --single --union -l 1 --force  --sinksource_method -c 1
python3 -u predict.py -o ../outfiles/SZ_1-layer_10-sinksource --single --union -l 1 --force  --sinksource_method -c 10
python3 -u predict.py -o ../outfiles/SZ_1-layer_50-sinksource --single --union -l 1 --force  --sinksource_method -c 50

python3 -u predict.py -o ../outfiles/ASD_1-layer_0-sinksource --single --union -l 1 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.01-sinksource --single --union -l 1 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_0.1-sinksource --single --union -l 1 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_1-sinksource --single --union -l 1 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_10-sinksource --single --union -l 1 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_1-layer_50-sinksource --single --union -l 1 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

## LAYER 2 SINGLE RUNS
python3 -u predict.py -o ../outfiles/SZ_2-layer_0-sinksource --single --union -l 2 --force  --sinksource_method -c 0
python3 -u predict.py -o ../outfiles/SZ_2-layer_0.01-sinksource --single --union -l 2 --force  --sinksource_method -c 0.01
python3 -u predict.py -o ../outfiles/SZ_2-layer_0.1-sinksource --single --union -l 2 --force  --sinksource_method -c 0.1
python3 -u predict.py -o ../outfiles/SZ_2-layer_1-sinksource --single --union -l 2 --force  --sinksource_method -c 1
python3 -u predict.py -o ../outfiles/SZ_2-layer_10-sinksource --single --union -l 2 --force  --sinksource_method -c 10
python3 -u predict.py -o ../outfiles/SZ_2-layer_50-sinksource --single --union -l 2 --force  --sinksource_method -c 50

python3 -u predict.py -o ../outfiles/ASD_2-layer_0-sinksource --single --union -l 2 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_0.01-sinksource --single --union -l 2 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_0.1-sinksource --single --union -l 2 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_1-sinksource --single --union -l 2 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_10-sinksource --single --union -l 2 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_2-layer_50-sinksource --single --union -l 2 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

## LAYER 3 SINGLE RUNS
python3 -u predict.py -o ../outfiles/SZ_3-layer_0-sinksource --single --union -l 3 --force  --sinksource_method -c 0
python3 -u predict.py -o ../outfiles/SZ_3-layer_0.01-sinksource --single --union -l 3 --force  --sinksource_method -c 0.01
python3 -u predict.py -o ../outfiles/SZ_3-layer_0.1-sinksource --single --union -l 3 --force  --sinksource_method -c 0.1
python3 -u predict.py -o ../outfiles/SZ_3-layer_1-sinksource --single --union -l 3 --force  --sinksource_method -c 1
python3 -u predict.py -o ../outfiles/SZ_3-layer_10-sinksource --single --union -l 3 --force  --sinksource_method -c 10
python3 -u predict.py -o ../outfiles/SZ_3-layer_50-sinksource --single --union -l 3 --force  --sinksource_method -c 50

python3 -u predict.py -o ../outfiles/ASD_3-layer_0-sinksource --single --union -l 3 --force  --sinksource_method -c 0 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_0.01-sinksource --single --union -l 3 --force  --sinksource_method -c 0.01 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_0.1-sinksource --single --union -l 3 --force  --sinksource_method -c 0.1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_1-sinksource --single --union -l 3 --force  --sinksource_method -c 1 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_10-sinksource --single --union -l 3 --force  --sinksource_method -c 10 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt
python3 -u predict.py -o ../outfiles/ASD_3-layer_50-sinksource --single --union -l 3 --force  --sinksource_method -c 50 -d ../infiles/ASD_positives.txt -n ../infiles/ASD_negatives.txt

echo DONE.