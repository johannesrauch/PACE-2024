#!/bin/sh
# installs the highs lib locally in highs folder
cd highs
rm -rf build
mkdir build
cd build
cmake ..
cmake --build .
cmake --install .
