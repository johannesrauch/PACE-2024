# PACE
This project includes code for [PACE 2024](https://pacechallenge.org/2024/).
That is, this project implements exact solvers and heuristics to solve the One-Sided Crossing Minimization problem.
Its main executable is the exact solver `weberknecht`.

## Dependencies
- Modern C++ compiler; we use `c++17`.
- `CMake`, version >= 3.15.
- [HiGHS](https://github.com/ERGO-Code/HiGHS).
- [`fmt/printf.hpp`](https://github.com/afborchert/fmt/blob/master/printf.hpp) (copy included in this project).

## Building
Currently only linux builds are supported.

### Installing HiGHS
This works automatically, just follow the steps below.

### Building `weberknecht`
The command-chain is
```
mkdir build
cd build
cmake ..
cmake --build .
```
This creates the executable `weberknecht` in `build`, which is the exact solver.

### Building Tests

For testing, execute
```
cmake -DPACE_BUILD_TESTS=ON ..
```
instead of the third command.

## Contributing
You are welcome to contribute.
Please follow the following rules.

### Formatting
Please use clang-format with
```
{ BasedOnStyle: Google, IndentWidth: 4, ColumnLimit: 120 }
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
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html).

## Project Status
Mostly done; some minor work in progress.

### Possible Future Work
Unfortunately, HiGHS does not (yet) support lazy constraint callbacks in its MIP solver.
Therefore, this project implements a specialized branch and cut MIP solver using the HiGHS LP solver, which is not ideal.
Once this feature is introduced in HiGHS, it might be worthwile to use it.
