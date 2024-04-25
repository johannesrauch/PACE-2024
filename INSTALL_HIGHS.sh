#!/bin/sh
# installs the highs lib locally in highs folder
cd highs
rm -rf build include lib
mkdir build
cd build
cmake ..
cmake --build .
cmake --install .
