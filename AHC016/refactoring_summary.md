# AHC016 Refactoring Summary

## Completed Refactorings

### 1. Array Type Modernization
- Replaced raw C-style arrays with `std::array` for fixed-size arrays
- Replaced raw arrays with `std::vector` for dynamic arrays
- Created template functions to handle bitset operations with different sizes

### 2. Function Name Refactoring (Minimal)
- `RandmizeGraph` → `ApplyNoiseToGraph`
- `InitB` → `ReceiveNoisyGraph`
- `ComputeAnswer` → `DecodeGraphIndex`
- Added Japanese comments explaining variable meanings

### 3. Common Pattern Extraction

#### MAX_ATTEMPTS Loop Pattern
- Extracted clique-finding loops into `findClique()` function
- Handles cliques of size 3, 4, or 5
- Replaced 9 occurrences, saving ~200 lines

#### Greedy Elimination Pattern
- Created `performGreedyElimination()` function
- Replaced 12 occurrences, saving 136 lines

#### Randomized Greedy Elimination
- Created `performRandomizedGreedyElimination()` function
- Replaced 4 occurrences, saving 96 lines

#### Solver Helper Functions
- Created 6 helper functions for common Solver patterns:
  - `initializeCountAndFlags()`
  - `findBestMatchingSingleCore()`
  - `findBestMatchingPairCores()`
  - `recountForSubset()`
  - `performSecondCoreSearch()`
  - `performBitsetOptimization()`
- Refactored Solver2, Solver4, Solver5, and Solver6

### 4. Magic Number Extraction
- Added named constants:
  - `ARRAY_SIZE = 1000`
  - `SPECIAL_MODE = 1000`
  - `INITIAL_DIFF = 1000`
  - `FLIP_ITERATIONS = 1000`
  - `TIME_LIMIT_MS = 360000.0`
  - `MISMATCH_PENALTY = 0.9`

### 5. Function Pointer Array
- Replaced 24-case if-else chain in `DecodeGraphIndex()` with function pointer array
- Cleaner, more maintainable code

## Total Impact
- Reduced code duplication by ~400+ lines
- Improved code maintainability and readability
- Maintained identical functionality throughout
- All refactorings tested and verified to produce same output

## Not Refactored
- InitNumArray functions (too diverse to consolidate effectively)
- solve() function (412 lines, but complex mode-dependent behavior)