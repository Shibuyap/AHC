#!/usr/bin/env python3

# Read the file
with open('AHC013_backup.cpp', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find key positions
innermethod_start = None
first_edge_creation = None
second_edge_creation = None
innermethod_end = None
direction_functions_start = None

for i, line in enumerate(lines):
    if 'int InnerMethod(' in line:
        innermethod_start = i
    if innermethod_start and '// 新たな辺の作成' in line:
        if first_edge_creation is None:
            first_edge_creation = i
        else:
            second_edge_creation = i
    if innermethod_start and line.strip() == 'return isDo;':
        # Find the closing brace after return
        for j in range(i+1, len(lines)):
            if lines[j].strip() == '}':
                innermethod_end = j
                break
    if 'void processRightMove(' in line:
        direction_functions_start = i

print(f"InnerMethod start: {innermethod_start}")
print(f"First edge creation: {first_edge_creation}")
print(f"Second edge creation: {second_edge_creation}")
print(f"InnerMethod end: {innermethod_end}")
print(f"Direction functions start: {direction_functions_start}")

# Build the new file
new_lines = []

# Copy everything up to the first edge creation setup
for i in range(0, 1149):
    new_lines.append(lines[i])

# Add the direction function calls and then jump to edge creation
new_lines.append("  // The direction-specific functions handle all the movement logic\n")
new_lines.append("  // Continue with edge creation and remaining logic\n")
new_lines.append("\n")

# Copy from the second edge creation to the end of InnerMethod
for i in range(second_edge_creation, innermethod_end + 1):
    new_lines.append(lines[i])

# Copy everything after InnerMethod
for i in range(innermethod_end + 1, len(lines)):
    new_lines.append(lines[i])

# Write the fixed file
with open('AHC013.cpp', 'w', encoding='utf-8') as f:
    f.writelines(new_lines)

print("Fixed AHC013.cpp")