wget http://www.martinwerner.de/files/dataset-sample.tgz
tar -xzvf dataset-sample.tgz
python3 scripts/giscup_to_mpk.py files/ SanFrancisco.mpk
rm dataset-sample.tgz
rm -r files/
python3 scripts/split_dataset.py SanFrancisco.mpk 100
