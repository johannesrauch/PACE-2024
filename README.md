# PACE2024
## Dependencies
- Modern version of `gcc`; we use `c++17`.
- `make`.
- [GLPK](https://www.gnu.org/software/glpk/). The makefile is configured such that the glpk library should be in the standard folder `/usr/local/include`.

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
for all tests (except the solver tests with large instances).
The command
```
make clean
```
deletes all executables in the test folder.

## Formatting
Please use clang-format with
```
{ BasedOnStyle: Google, IndentWidth: 4, ColumnLimit: 0 }
```
options.

Please use 
```
#ifndef PACE2024_[YOUR FILE NAME IN CAPS WITH UNDERSCORES]_HPP
#define PACE2024_[YOUR FILE NAME IN CAPS WITH UNDERSCORES]_HPP

// your code here

#endif
```
to avoid multiple header inclusions.
Please use the namespace `pace2024` for the project.

## Add files
```
cd existing_repo
git remote add origin https://gitlab.uni-ulm.de/zhh72/pace2024.git
git branch -M main
git push -uf origin main
```

## Roadmap
- [ ] Switch linear program solver to [CLP](https://github.com/coin-or/Clp).
- [ ] Implement improvements (todos) mentioned in `branch_and_cut::solve`.

## Authors
Johannes Rauch.

## License
[MIT License](https://en.wikipedia.org/wiki/MIT_License).

## Project status
Work in progress.
