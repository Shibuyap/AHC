#!/usr/bin/env python3
import re
import sys

def replace_rep_macros(file_path):
    with open(file_path, 'r', encoding='utf-8-sig') as f:
        content = f.read()
    
    # Pattern to match rep(variable, limit)
    rep_pattern = r'\brep\s*\(\s*([^,]+?)\s*,\s*([^)]+?)\s*\)'
    # Pattern to match srep(variable, start, limit)
    srep_pattern = r'\bsrep\s*\(\s*([^,]+?)\s*,\s*([^,]+?)\s*,\s*([^)]+?)\s*\)'
    # Pattern to match drep(variable, limit)
    drep_pattern = r'\bdrep\s*\(\s*([^,]+?)\s*,\s*([^)]+?)\s*\)'
    
    # Replace rep(i, n) with for (int i = 0; i < (n); ++i)
    content = re.sub(rep_pattern, r'for (int \1 = 0; \1 < (\2); ++\1)', content)
    
    # Replace srep(i, s, t) with for (int i = s; i < t; ++i)
    content = re.sub(srep_pattern, r'for (int \1 = \2; \1 < \3; ++\1)', content)
    
    # Replace drep(i, n) with for (int i = (n)-1; i >= 0; --i)
    content = re.sub(drep_pattern, r'for (int \1 = (\2)-1; \1 >= 0; --\1)', content)
    
    # Write the result back with UTF-8 BOM
    with open(file_path, 'w', encoding='utf-8-sig') as f:
        f.write(content)
    
    print(f"Replaced macros in {file_path}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        replace_rep_macros(sys.argv[1])
    else:
        print("Usage: python replace_rep_macros.py <file_path>")