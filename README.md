# frechet_simsearch
A C++17 implementation of data structures for approximate trajectory similarity search under Fréchet distance, used in the experiments of the paper, Shunsuke Kanda, Koh Takeuchi, Keisuke Fujii, and Yasuo Tabei "[Succinct Trit-array Trie for Scalable Trajectory Similarity Search](https://arxiv.org/abs/2005.10917)," In *28th ACM SIGSPATIAL*, 2020.

## Build instructions

You can download and compile this library by the following commands:

```shell
$ git clone https://github.com/kampersanda/frechet_simsearch.git
$ cd frechet_simsearch
$ mkdir build
$ cd build
$ cmake ..
$ make -j
```

After the commands, the executables will be produced in `build/bin` directory.

The code is written in C++17, so please install g++ >= 7.0 or clang >= 4.0. The following dependencies have to be installed to compile the library: CMake >= 3.0 (for the build system), Boost >= 1.42, [sdsl-lite](https://github.com/simongog/sdsl-lite), and [msgpack-c](https://github.com/msgpack/msgpack-c).

## Examples to produce trajectory datasets

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

Then, data set `SanFrancisco_base.mpk` and query set `SanFrancisco_query.mpk` will be generated in the msgpack format.

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

Then, data set `PortoTaxi_base.mpk` and query set `PortoTaxi_query.mpk` will be generated in the msgpack format.

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

An example to analyze tSTAT on `PortoTaxi` using Fréchet range 7000 is as follows.

The LSH performance (using default parameters) can be analyzed in the following commands:

```shell
$ ./bin/make_range_groundtruth -b PortoTaxi_base.mpk -q PortoTaxi_query.mpk -r 7000
$ ./bin/analyze_array_score -b PortoTaxi_base.mpk -q PortoTaxi_query.mpk -r 7000
```

The search performance (using default parameters) can be analyzed in the following commands:

```shell
$ ./bin/analyze_mtrie_perf -b PortoTaxi_base.mpk -q PortoTaxi_query.mpk -r 7000
```

## Licensing

This program is available for only academic use, basically. For the academic use, please keep MIT License. For the commercial use, please keep GPL 2.0 and make a contact to [Shunsuke Kanda](shnsk.knd@gmail.com).

If you use the library, please cite the following paper:

```tex
@inproceedings{kanda2020tstat,
  author = {Kanda, Shunsuke and Takeuchi, Koh and Fujii, Keisuke and Tabei, Yasuo},
  title = {Succinct trit-array trie for scalable trajectory similarity search},
  booktitle = {Proceedings of the 28th ACM SIGSPATIAL International Conference on Advances in Geographic Information Systems (SIGSPATIAL)},
  year = {2020}
}
```

