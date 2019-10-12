# frechet_simsearch
[![experimental](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)

This is an experimental library for trajectory similarity search under Fréchet distance.

## Build instructions

```shell
$ git clone https://github.com/kampersanda/frechet_simsearch.git
$ cd frechet_simsearch
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## External libraries

You have to install the external libraries:

- [Boost](https://www.boost.org)
- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [msgpack-c](https://github.com/msgpack/msgpack-c)

## Examples to make trajectory datasets

### [GISCUP 2017 dataset](http://sigspatial2017.sigspatial.org/giscup2017/download)

```shell
$ wget http://www.martinwerner.de/files/dataset-sample.tgz
$ tar -xzvf dataset-sample.tgz
$ mkdir data
$ python3 scripts/giscup_to_mpk.py files/ data/SanFrancisco.mpk
$ rm dataset-sample.tgz
$ rm -r files/
$ python3 scripts/split_dataset.py SanFrancisco.mpk 100
```

The statistics are

```shell
$ ./scripts/analyze_curve_statistics.py SanFrancisco.mpk 
num_curves: 20199
num_points: 5007230
min_length: 10
max_length: 768
mean_length: 247.8949452943215
median_length: 222.0
```

### [ECML/PKDD 2015 dataset](https://archive.ics.uci.edu/ml/datasets/Taxi+Service+Trajectory+-+Prediction+Challenge,+ECML+PKDD+2015#)

```shell
$ wget https://archive.ics.uci.edu/ml/machine-learning-databases/00339/train.csv.zip
$ unzip train.csv.zip
$ python3 scripts/porto_to_mpk.py train.csv PortoTaxi.mpk
$ rm train.csv
$ rm train.csv.zip
$ python3 scripts/split_dataset.py PortoTaxi.mpk 1000
```

The statistics are

```shell
$ ./scripts/analyze_curve_statistics.py PortoTaxi.mpk 
num_curves: 1704769
num_points: 83409386
min_length: 1
max_length: 3881
mean_length: 48.927089828592614
median_length: 41.0
```

## Licensing

This library is free software provided under [Apache License 2.0](https://github.com/kampersanda/frechet_simsearch/blob/master/LICENSE), following the License of [Cecca/FRESH](https://github.com/Cecca/FRESH). Our modifications are put in each source file.