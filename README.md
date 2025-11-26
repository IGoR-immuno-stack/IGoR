# IGoR: Inference and Generation Of Repertoires

This project is a fork of the original one available on github: https://github.com/qmarcou/IGoR.

The goal is this repository is to present what could be a new organization of the source code and the build tree in order to make it easier to use and modify.

## Prerequisites

Either `conda` or `mamba` (https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) is required in order to use the code properly. 

`Mamba` is just a c++ reimplementation of `conda` which provides better performances.

## Compilation

### MacOS

```shell
$ cd path/to/igor
$ mkdir build; cd build
$ mamba env update -f ../pkg/env/igor-macOS.yaml
$ mamba activate igor
$ cmake .. -DCMAKE_CXX_COMPILER=${CONDA_PREFIX}/bin/clang++ \
           -DCMAKE_C_COMPILER=${CONDA_PREFIX}/bin/clang
$ make
```

### Linux

```shell
$ cd path/to/igor
$ mkdir build; cd build
$ mamba env update -f ../pkg/env/igor.yaml
$ mamba activate igor
$ cmake .. -DCMAKE_CXX_COMPILER=${CONDA_PREFIX}/bin/g++ \
           -DCMAKE_C_COMPILER=${CONDA_PREFIX}/bin/gcc
$ make
```

## Example of Use

Firstly, it is necessary to clone `igor-models` repository to execute the tests. 

```shell
$ cd path/to/igor/..
$ git clone git@gitlab.inria.fr:ant-men/igor-stack/igor-models.git
```

Then, one can execute igor as follows:
```shell
$ cd path/to/igor
$ cd script
$ ./igor-compute_pgen human beta actcagctttgtatttctgtgccagcagcgtagattgggacagggggcctcctacgagcagtacgtcgggccg

```

```shell
tr: Illegal byte sequence
human beta /var/folders/zq/84q829qd5lbb_397bmzy51g40000gn/T/tmp.saaBn1Sz /var/folders/zq/84q829qd5lbb_397bmzy51g40000gn/T/tmp.saaBn1Sz/igorSeq.XXXXXXXXXXXXXXXXXXXXXX.txt
Batch name set to: #Q_
TXT extension detected for the input sequence file
Working directory set to: "/var/folders/zq/84q829qd5lbb_397bmzy51g40000gn/T/tmp.saaBn1Sz/"
Batch name set to: #Q_
Species parameter set to: human
Chain parameter set to: beta
Working directory set to: "/var/folders/zq/84q829qd5lbb_397bmzy51g40000gn/T/tmp.saaBn1Sz/"
Performing V alignments....
V_gene alignments [||||||||||||||||||||||||||||||||||||||||||||||||||]  Done.
Performing D alignments....
D_gene alignments [||||||||||||||||||||||||||||||||||||||||||||||||||]  Done.
Performing J alignments....
J_gene alignments [||||||||||||||||||||||||||||||||||||||||||||||||||]  Done.
Performing CDR3 sequence extraction ....
Batch name set to: #Q_
Species parameter set to: human
Chain parameter set to: beta
Working directory set to: "/var/folders/zq/84q829qd5lbb_397bmzy51g40000gn/T/tmp.saaBn1Sz/"
Read some model parms
GeneChoice read
GeneChoice read
GeneChoice read
Deletion read
Deletion read
Deletion read
Deletion read
Insertion read
Insertion read
DinucMarkov read
DinucMarkov read
Performing Evaluate/Inference iteration 1
Initializing probability bounds...
Initialization of probability bounds over.
Iteration 1 [||||||||||||||||||||||||||||||||||||||||||||||||||]  Done.
8.38274e-12
```

