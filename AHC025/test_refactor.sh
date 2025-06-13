#!/bin/bash

# Compile both versions
echo "Compiling original code..."
g++ -O2 -std=c++17 -o AHC025_original AHC025_backup.cpp
echo "Compiling refactored code..."
g++ -O2 -std=c++17 -o AHC025_refactored AHC025.cpp

# Create directories for output
mkdir -p original_output refactored_output

# Counter for matching results
matches=0

# Test on 100 cases
for i in {0..99}; do
    # Format case number with leading zeros
    case_num=$(printf "%04d" $i)
    
    # Run original
    ./AHC025_original < in/${case_num}.txt > original_output/${case_num}.txt 2>&1
    original_score=$(tail -n 1 original_output/${case_num}.txt | grep -o '[0-9]*')
    
    # Run refactored
    ./AHC025_refactored < in/${case_num}.txt > refactored_output/${case_num}.txt 2>&1
    refactored_score=$(tail -n 1 refactored_output/${case_num}.txt | grep -o '[0-9]*')
    
    # Compare scores
    if [ "$original_score" = "$refactored_score" ]; then
        ((matches++))
        echo "Case $case_num: MATCH (score: $original_score)"
    else
        echo "Case $case_num: DIFFER (original: $original_score, refactored: $refactored_score)"
    fi
done

echo ""
echo "Summary: $matches out of 100 cases matched"
if [ $matches -ge 80 ]; then
    echo "SUCCESS: At least 80% of cases match!"
else
    echo "FAILURE: Less than 80% of cases match"
fi