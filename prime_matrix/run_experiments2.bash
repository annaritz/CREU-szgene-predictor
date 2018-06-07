1set -e  ## this will mean that the script continues

echo SINGLE RUNS
python3 predict2.py --single --force --sinksource_constant 10 -l 1
python3 predict2.py --single --force --sinksource_constant 1 -l 2
python3 predict2.py --single --force --sinksource_constant 1 -l 3
python3 predict2.py --single --force --sinksource_constant 1 -l 4








