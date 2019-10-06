wget https://archive.ics.uci.edu/ml/machine-learning-databases/00339/train.csv.zip
unzip train.csv.zip
python3 scripts/porto_to_mpk.py train.csv PortoTaxi.mpk
rm train.csv
rm train.csv.zip
python3 scripts/split_dataset.py PortoTaxi.mpk 1000
