# PACE 2024
This repository provides exact and heuristic solvers for the One-Sided Crossing Minimization problem of [PACE 2024](https://pacechallenge.org/2024/).

Its main executables are the exact solver `weberknecht` (for the exact and parameterized track) and the heuristic solver `weberknecht_h` (for the heuristic track).

The solver description is `docs/main.pdf`.

## Dependencies
- Modern C++ compiler, we use the `c++17` language standard.
- `cmake`, version >= 3.15.
- [HiGHS](https://github.com/ERGO-Code/HiGHS) (MIT License).
- [`fmt/printf.hpp`](https://github.com/afborchert/fmt/blob/master/printf.hpp) (MIT License, copy included in this project).

## Building
Currently only linux builds are supported.
Windows and macOS builds may work too, but we did not test it.

### Installing HiGHS
`cmake` installs [HiGHS](https://github.com/ERGO-Code/HiGHS) automatically in the project folder.
Note that `cmake ..` may take a while since it clones the repository.

### Building the executables `weberknecht` and `weberknecht_h`
The command-chain is
```
mkdir build
cd build
cmake ..
cmake --build .
```
This creates the executables `weberknecht` and `weberknecht_h` in `build`.

### Building the tests

For testing, execute
```
cmake -DPACE_BUILD_TESTS=ON ..
```
instead `cmake ..`.
The test executables are built in `build/test`.

### Building with debug print

For testing, execute
```
cmake -DPACE_DEBUG_PRINT=ON ..
```
instead `cmake ..`.

## Running

After building, the executables are in the `build` folder.
Both `weberknecht` and `weberknecht_h` expect input over stdio, e.g.
```
./weberknecht < path/to/instance/42.gr
```

## Contributing
You are welcome to contribute.
Please follow the following rules.

### Formatting
Please use clang-format with
```
{ BasedOnStyle: Google, IndentWidth: 4, ColumnLimit: 80 }
```
options.

### Guard clauses
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
See `LICENSE.md`.

## Project Status
Mostly done; some minor work in progress.

Some workable points:
- Unfortunately, HiGHS does not (yet) support lazy constraint callbacks in its MIP solver.
Therefore, this project implements a specialized branch and cut MIP solver using the HiGHS LP solver, which is not ideal.
Once this feature is introduced in HiGHS, it might be worthwile to use it.
