#!/usr/bin/env python3

# Read the file
with open('AHC013.cpp', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find the lines to remove
# We want to keep up to and including the first "// 新たな辺の作成" section setup
# Then skip to the actual edge creation logic at the second occurrence

new_lines = []
i = 0
while i < len(lines):
    if i < 1166:  # Keep everything before the old direction code
        new_lines.append(lines[i])
    elif i >= 2100:  # Keep everything from the second "// 新たな辺の作成" onward
        new_lines.append(lines[i])
    i += 1

# Write the fixed file
with open('AHC013_fixed.cpp', 'w', encoding='utf-8') as f:
    f.writelines(new_lines)

print("Created AHC013_fixed.cpp with the corrected InnerMethod function")