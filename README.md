# PACE
## Dependencies
- Modern version of `gcc`; we use `c++17`.
- `make`.
- [HiGHS](https://github.com/ERGO-Code/HiGHS). `make` expects the headers in `/usr/local/include/highs` and the library in `/usr/local/lib`.

## Compiling
### Tests
The command
```
make test.bipartite_graph
```
runs the bipartite_graph.hpp test, for instance.
Run 
```
make test.all
```
for all tests (on `tiny_test_set`).
The command
```
make clean
```
deletes all executables and logs in the test folder.

## Formatting
Please use clang-format with
```
{ BasedOnStyle: Google, IndentWidth: 4, ColumnLimit: 120 }
```
options.

Please use 
```
#ifndef PACE_[FILEPATH IN CAPS WITH UNDERSCORES]_HPP
#define PACE_[FILEPATH IN CAPS WITH UNDERSCORES]_HPP

// your code here

#endif
```
to avoid multiple header inclusions.
Please use the namespace `pace` for the project.

## Authors
Johannes Rauch.

## License
[MIT License](https://en.wikipedia.org/wiki/MIT_License).

## Project status
Work in progress.
