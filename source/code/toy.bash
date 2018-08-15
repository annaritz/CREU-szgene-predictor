## no SS
python3 predict.py --toy -l 1 --force
python3 predict.py --toy -l 2 --force

## SS
python3 predict.py --toy -l 1 --force --sinksource_method -c 0
python3 predict.py --toy -l 1 --force --sinksource_method -c 0.01
python3 predict.py --toy -l 1 --force --sinksource_method -c 0.1
python3 predict.py --toy -l 1 --force --sinksource_method -c 1
python3 predict.py --toy -l 1 --force --sinksource_method -c 10

python3 predict.py --toy -l 2 --force --sinksource_method -c 0
python3 predict.py --toy -l 2 --force --sinksource_method -c 0.01
python3 predict.py --toy -l 2 --force --sinksource_method -c 0.1
python3 predict.py --toy -l 2 --force --sinksource_method -c 1
python3 predict.py --toy -l 2 --force --sinksource_method -c 10