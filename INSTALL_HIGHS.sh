#!/bin/sh
cd highs
rm -rf build
mkdir build
cd build
cmake ..
cmake --build .
cmake --install .
