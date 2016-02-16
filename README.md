Instruction to set up fourtops can be found here:

https://docs.google.com/document/d/1cbFPpLJSKc3ketyKBFjR1N0IuiJh-PmRvIorqoeVBqU/edit


This package utilises glog library for logging purposes. Therefore pull glog from github and install it.
So far only v0.3.3 had no installation problems at m-machines

```sh
git clone https://github.com/google/glog.git
git checkout tags/v0.3.3
./configure --prefix=$PWD
make
make install
```

To compile SingleLepton lib
```sh
mkdir build
cd build
cmake -DTopBrussels_SOURCE_DIR=/user/direcrory/TopBrussels ..
make -j
make install
```
