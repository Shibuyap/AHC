#!/bin/bash

# Test on 20 cases
matches=0
total=20

echo "Testing refactored code on $total cases..."

for i in $(seq 0 19); do
    case_num=$(printf "%04d" $i)
    
    # Run original
    original_score=$(./AHC025_original < in/${case_num}.txt 2>&1 | tail -1 | grep -oP 'score =\s*\K[0-9]+')
    
    # Run refactored
    refactored_score=$(./AHC025_refactored < in/${case_num}.txt 2>&1 | tail -1 | grep -oP 'score =\s*\K[0-9]+')
    
    # Compare scores
    if [ "$original_score" = "$refactored_score" ]; then
        ((matches++))
        echo "Case $case_num: MATCH (score: $original_score)"
    else
        echo "Case $case_num: DIFFER (original: $original_score, refactored: $refactored_score)"
    fi
done

echo ""
echo "Summary: $matches out of $total cases matched"
percent=$((matches * 100 / total))
echo "Match rate: $percent%"