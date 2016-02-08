Instruction to set up fourtops can be found here:

https://docs.google.com/document/d/1cbFPpLJSKc3ketyKBFjR1N0IuiJh-PmRvIorqoeVBqU/edit


To compile SingleLepton lib

mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=~/lib ..
make -j
make install

