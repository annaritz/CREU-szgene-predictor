set -e  ## this will mean that the script continues

echo SINGLE RUNS
#python3 predict.py -g networkfiles/brain_top_geq_0.500.txt -o outfiles/network_0.500 --matrix --plot --roc --format
#python3 predict.py -g networkfiles/brain_top_geq_0.400.txt -o outfiles/network_0.400 --matrix --plot --roc --format
#python3 predict.py -g networkfiles/brain_top_geq_0.300.txt -o outfiles/network_0.300 --matrix --plot --roc --format
#python3 predict.py -g networkfiles/brain_top_geq_0.200.txt -o outfiles/network_0.200 --matrix --plot --roc --format
#python3 predict.py -g networkfiles/brain_top_geq_0.150.txt -o outfiles/network_0.150 --matrix --plot --roc --format
#python3 predict.py -g networkfiles/brain_top_geq_0.125.txt -o outfiles/network_0.125 --matrix --plot --roc --format
#python3 predict.py -g networkfiles/brain_top_geq_0.100.txt -o outfiles/network_0.100 --matrix --plot --roc --format

echo AUCs
#python3 predict.py -g networkfiles/brain_top_geq_0.500.txt -o outfiles/network_0.500 --matrix --auc
python3 predict.py -g networkfiles/brain_top_geq_0.400.txt -o outfiles/network_0.400 --matrix --auc
python3 predict.py -g networkfiles/brain_top_geq_0.300.txt -o outfiles/network_0.300 --matrix --auc
python3 predict.py -g networkfiles/brain_top_geq_0.200.txt -o outfiles/network_0.200 --matrix --auc
python3 predict.py -g networkfiles/brain_top_geq_0.150.txt -o outfiles/network_0.150 --matrix --auc
python3 predict.py -g networkfiles/brain_top_geq_0.125.txt -o outfiles/network_0.125 --matrix --auc
python3 predict.py -g networkfiles/brain_top_geq_0.100.txt -o outfiles/network_0.100 --matrix --auc

## aggregate
python3 predict.py --aggregate outfiles/network_0.125 --aggregate outfiles/network_0.150 --aggregate outfiles/network_0.200 --aggregate outfiles/network_0.300 --aggregate outfiles/network_0.400 --aggregate outfiles/network_0.500 --aggregate_names 0.125 --aggregate_names 0.150 --aggregate_names 0.200 --aggregate_names 0.300 --aggregate_names 0.400 --aggregate_names 0.500
echo DONE.
