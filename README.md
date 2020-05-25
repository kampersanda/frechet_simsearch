# frechet_simsearch
[![experimental](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)

This is an experimental library for approximate trajectory similarity search under Fréchet distance, used in the paper [Succinct Trit-array Trie for Scalable Trajectory Similarity Search](https://arxiv.org/abs/2005.10917).

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

## Example to make trajectory datasets

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

Then, data set `SanFrancisco_base.mpk` and query set `SanFrancisco_query.mpk` will be generated.

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

Then, data set `PortoTaxi_base.mpk` and query set `PortoTaxi_query.mpk` will be generated.

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

## Example to analyze search methods

An example to analyze tSTAT on `PortoTaxi` using Fréchet range 7263 is as follows.

The LSH performance (using other default parameters) can be analyzed in the following commands:

```shell
$ ./bin/make_range_groundtruth -b PortoTaxi_base.mpk -q PortoTaxi_query.mpk -r 7263
$ ./bin/analyze_array_score -b PortoTaxi_base.mpk -q PortoTaxi_query.mpk -r 7263
```

The search performance (using other default parameters) can be analyzed in the following commands:

```shell
$ ./bin/analyze_mtrie_perf -b PortoTaxi_base.mpk -q PortoTaxi_query.mpk -r 7263
```

## Licensing

This library is free software provided under [Apache License 2.0](https://github.com/kampersanda/frechet_simsearch/blob/master/LICENSE), following the License of [Cecca/FRESH](https://github.com/Cecca/FRESH). Our modifications are put in each source file.