# Clique Loop Refactoring Script

## Overview
The `refactor_clique_loops.py` script automatically refactors repetitive clique-finding loops in AHC016.cpp by replacing them with calls to the `findClique` function.

## What it does

1. **Identifies MAX_ATTEMPTS loops**: Finds all loops that follow the pattern:
   ```cpp
   for (int loop1 = 0; loop1 < MAX_ATTEMPTS; ++loop1) {
     // clique finding code
   }
   ```

2. **Analyzes each loop** to determine:
   - Clique size (3, 4, or 5)
   - Mark value (1 or 2)
   - Which cores vector is used (cores, cores1, or cores2)

3. **Replaces loops** with appropriate `findClique` calls:
   ```cpp
   if (findClique(kouho, f, cores2, 5, 2)) {
     // code that was after the loop
   }
   ```

4. **Handles conditional logic**: If the original code had an `if (cores.size() > 0)` check after the loop, it merges that logic into the replacement.

## Results

Running on AHC016.cpp:
- Found and replaced 10 clique-finding loops
- Reduced file from 4305 lines to 4163 lines (saved 142 lines)
- Replaced patterns:
  - 5-clique loops with mark=1 and mark=2
  - 4-clique loops with mark=1
  - 3-clique loops with mark=2

## Usage

```bash
python3 refactor_clique_loops.py
```

The script reads `AHC016.cpp` and creates `AHC016_refactored.cpp` with the refactored code.

## Technical Details

The script uses regular expressions to:
1. Find loop patterns
2. Extract loop bodies
3. Analyze clique parameters
4. Generate appropriate replacements

It handles edge cases like:
- Different variable names for candidate arrays (kouho)
- Various cores vector names (cores, cores1, cores2)
- Conditional logic after loops
- Proper indentation preservation