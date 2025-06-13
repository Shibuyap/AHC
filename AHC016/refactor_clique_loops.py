#!/usr/bin/env python3
"""
Script to refactor AHC016.cpp by replacing MAX_ATTEMPTS loops with findClique function calls.
"""

import re
import sys
from pathlib import Path


def extract_clique_loop(content, start_pos):
    """Extract the entire clique-finding loop starting at start_pos."""
    # Find the matching closing brace
    brace_count = 0
    i = start_pos
    loop_start = start_pos
    
    # Find the opening brace of the loop
    while i < len(content) and content[i] != '{':
        i += 1
    if i >= len(content):
        return None, None
    
    i += 1  # Skip the opening brace
    brace_count = 1
    
    while i < len(content) and brace_count > 0:
        if content[i] == '{':
            brace_count += 1
        elif content[i] == '}':
            brace_count -= 1
        i += 1
    
    if brace_count != 0:
        return None, None
    
    loop_end = i
    return content[loop_start:loop_end], loop_end


def analyze_clique_loop(loop_content):
    """Analyze the loop to determine clique size and mark value."""
    # Look for the core array declaration
    core_match = re.search(r'(?:int|vector<int>)\s+core(?:\[(\d+)\]|\((\d+)\))', loop_content)
    if not core_match:
        return None, None, None
    
    clique_size = core_match.group(1) or core_match.group(2)
    if not clique_size:
        # Try to find it from the loop bounds
        size_match = re.search(r'for\s*\(\s*int\s+i\s*=\s*0\s*;\s*i\s*<\s*\(?\s*(\d+)\s*\)?', loop_content)
        if size_match:
            clique_size = size_match.group(1)
        else:
            return None, None, None
    
    # Find the mark value being assigned
    mark_match = re.search(r'f\[core\[i\]\]\s*=\s*(\d+)', loop_content)
    if not mark_match:
        return None, None, None
    
    mark_value = mark_match.group(1)
    
    # Check if this is followed by cores.push_back or cores1/cores2.push_back
    cores_var = None
    if re.search(r'cores\.push_back', loop_content):
        cores_var = 'cores'
    elif re.search(r'cores1\.push_back', loop_content):
        cores_var = 'cores1'
    elif re.search(r'cores2\.push_back', loop_content):
        cores_var = 'cores2'
    
    return int(clique_size), int(mark_value), cores_var


def find_kouho_variable(content, loop_start):
    """Find the kouho variable used in the loop."""
    # Look backwards for kouho variable declaration/usage
    search_start = max(0, loop_start - 1000)
    search_content = content[search_start:loop_start]
    
    # Common patterns for kouho usage
    kouho_patterns = [
        r'vector<int>\s+(\w+kouho\w*)',
        r'(\w+kouho\w*)\[Rand',
        r'(\w+kouho\w*)\.size\(\)'
    ]
    
    kouho_var = None
    for pattern in kouho_patterns:
        match = re.search(pattern, search_content[::-1])
        if match:
            kouho_var = match.group(1)[::-1]
            break
    
    if not kouho_var:
        # Default to 'kouho'
        kouho_var = 'kouho'
    
    return kouho_var


def find_if_condition_after_loop(content, loop_end):
    """Find if there's an if-condition checking cores size after the loop."""
    # Look for pattern like: if (cores.size() > 0) or if (cores1.size() > 0)
    search_end = min(len(content), loop_end + 500)
    search_content = content[loop_end:search_end]
    
    if_match = re.search(r'if\s*\(\s*(cores\d?)\.size\(\)\s*>\s*0\s*\)', search_content)
    if if_match:
        return if_match.group(1)
    
    return None


def replace_clique_loops(content):
    """Replace all MAX_ATTEMPTS loops with findClique calls."""
    # Find all MAX_ATTEMPTS loops
    pattern = r'for\s*\(\s*int\s+loop1\s*=\s*0\s*;\s*loop1\s*<\s*MAX_ATTEMPTS\s*;\s*\+\+loop1\s*\)'
    
    replacements = []
    
    for match in re.finditer(pattern, content):
        start_pos = match.start()
        loop_content, loop_end = extract_clique_loop(content, start_pos)
        
        if not loop_content:
            continue
        
        clique_size, mark_value, cores_var = analyze_clique_loop(loop_content)
        
        if clique_size is None or mark_value is None:
            print(f"Warning: Could not analyze loop at position {start_pos}")
            continue
        
        print(f"Found {clique_size}-clique loop with mark={mark_value} at position {start_pos}")
        
        kouho_var = find_kouho_variable(content, start_pos)
        
        # Check if there's an if-condition after the loop
        cores_check = find_if_condition_after_loop(content, loop_end)
        
        # Build the replacement
        if cores_var:
            replacement = f"if (findClique({kouho_var}, f, {cores_var}, {clique_size}, {mark_value})) {{\n"
            
            # Find the code after the loop that should be included in the if block
            if cores_check:
                # Find the matching if block
                if_start = content.find(f"if ({cores_check}.size() > 0)", loop_end)
                if if_start != -1:
                    if_content, if_end = extract_clique_loop(content, if_start)
                    if if_content:
                        # Extract just the body of the if statement
                        body_match = re.search(r'\{(.+)\}', if_content, re.DOTALL)
                        if body_match:
                            replacement += body_match.group(1)
                        replacements.append((loop_end, if_end, ""))  # Remove the original if block
            
            replacement += "      }"
        else:
            # Simple replacement without cores variable - create a temporary one
            replacement = f"{{\n        vector<int> tempCores;\n        findClique({kouho_var}, f, tempCores, {clique_size}, {mark_value});\n      }}"
        
        replacements.append((start_pos, loop_end, replacement))
    
    # Apply replacements in reverse order to maintain positions
    replacements.sort(key=lambda x: x[0], reverse=True)
    
    result = content
    for start, end, replacement in replacements:
        result = result[:start] + replacement + result[end:]
    
    return result, len(replacements)


def main():
    input_file = Path("/mnt/c/Programming/AHC/AHC016/AHC016.cpp")
    output_file = Path("/mnt/c/Programming/AHC/AHC016/AHC016_refactored.cpp")
    
    if not input_file.exists():
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)
    
    print(f"Reading {input_file}...")
    content = input_file.read_text()
    
    print("Analyzing and replacing clique-finding loops...")
    refactored_content, num_replacements = replace_clique_loops(content)
    
    if num_replacements > 0:
        print(f"Writing refactored code to {output_file}...")
        output_file.write_text(refactored_content)
        print(f"Successfully replaced {num_replacements} clique-finding loops")
        
        # Show a summary of changes
        original_lines = content.count('\n')
        refactored_lines = refactored_content.count('\n')
        print(f"Original file: {original_lines} lines")
        print(f"Refactored file: {refactored_lines} lines")
        print(f"Lines saved: {original_lines - refactored_lines}")
    else:
        print("No clique-finding loops found to replace")


if __name__ == "__main__":
    main()