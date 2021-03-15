#!/bin/bash

source_dir=`pwd`
external_dir=${source_dir}/external
mkdir -p external
cd ${external_dir}
# build SZ (to use ZSTD compressor)
git clone https://github.com/szcompressor/SZ.git
cd SZ
git reset --hard f48d2f27a5470a28e900db9b46bb3344a2bc211f
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${external_dir}/SZ/install ..
make -j 8
make install

cd ${source_dir}
mkdir -p build
cd build
cmake ..
make -j 8
