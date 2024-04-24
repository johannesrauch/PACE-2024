#!/bin/sh
rm -rf highs/bin highs/build highs/include highs/lib
cd highs
mkdir build
cd build
cmake ..
rm -rf _deps/highs-src/.git
cd ../..
tar -cvz -f sub.tar.gz highs/build src CMakeLists.txt
