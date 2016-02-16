# Installation instructions

Instruction to set up fourtops can be found here:

https://docs.google.com/document/d/1cbFPpLJSKc3ketyKBFjR1N0IuiJh-PmRvIorqoeVBqU/edit


## glog library installation
This package utilises `glog` library for logging purposes. Therefore pull `glog` from the github and install it.
So far only v0.3.3 had no installation problems at m-machines

```sh
git clone https://github.com/google/glog.git
git checkout tags/v0.3.3
./configure --prefix=$PWD
make -j
make install
```

## gtest library installation
For testing (optional) google test framework is used. In order to install `googletest` do:

```sh
git clone https://github.com/google/googletest.git
cd googletest
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/googletest  ..
make -j
make install
```

## SingleLepton library and analysis code build instructions
This step assumes that a folder named `TopBrussels` containing `TopTreeAnalysisBase` and `TopTreeProducer` directories exists
To compile `SingleLepton` lib
```sh
mkdir build
cd build
cmake -DTopBrussels_SOURCE_DIR=/user/direcrory/TopBrussels ..
make -j
make install
```
