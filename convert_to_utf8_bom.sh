#!/bin/bash

# Convert all C++, header, and text files to UTF-8 with BOM

echo "Converting files to UTF-8 with BOM..."

# Function to add UTF-8 BOM if not present
add_bom_if_needed() {
    local file="$1"
    # Check if file starts with BOM
    if ! head -c 3 "$file" | hexdump -e '3/1 "%02x"' | grep -q "efbbbf"; then
        # File doesn't have BOM, add it
        echo "Adding BOM to: $file"
        # Create temp file with BOM
        printf '\xEF\xBB\xBF' > "${file}.tmp"
        cat "$file" >> "${file}.tmp"
        mv "${file}.tmp" "$file"
    fi
}

# Function to convert encoding and add BOM
convert_file() {
    local file="$1"
    local encoding=$(file -bi "$file" | sed 's/.*charset=//')
    
    echo "Processing: $file (detected encoding: $encoding)"
    
    # Skip binary files
    if file -b "$file" | grep -q "binary"; then
        echo "  Skipping binary file"
        return
    fi
    
    # Try to detect if it's Shift-JIS
    if [[ "$encoding" == "unknown-8bit" ]] || [[ "$encoding" == "iso-8859-1" ]]; then
        # Try converting from Shift-JIS
        if iconv -f SHIFT-JIS -t UTF-8 "$file" > "${file}.utf8.tmp" 2>/dev/null; then
            echo "  Converted from Shift-JIS to UTF-8"
            mv "${file}.utf8.tmp" "$file"
        else
            # Try converting from Windows-1252
            if iconv -f WINDOWS-1252 -t UTF-8 "$file" > "${file}.utf8.tmp" 2>/dev/null; then
                echo "  Converted from Windows-1252 to UTF-8"
                mv "${file}.utf8.tmp" "$file"
            else
                echo "  Warning: Could not convert encoding, keeping original"
                rm -f "${file}.utf8.tmp"
            fi
        fi
    elif [[ "$encoding" != "utf-8" ]] && [[ "$encoding" != "us-ascii" ]]; then
        # Try generic conversion
        if iconv -f "$encoding" -t UTF-8 "$file" > "${file}.utf8.tmp" 2>/dev/null; then
            echo "  Converted from $encoding to UTF-8"
            mv "${file}.utf8.tmp" "$file"
        else
            echo "  Warning: Could not convert from $encoding"
        fi
    fi
    
    # Add BOM if needed
    add_bom_if_needed "$file"
}

# Find and process all relevant files
find /mnt/c/Programming/AHC -type f \( -name "*.cpp" -o -name "*.h" -o -name "*.md" -o -name "*.txt" \) \
    -not -path "*/x64/*" -not -path "*/.vs/*" -not -path "*/Debug/*" -not -path "*/Release/*" | \
while read -r file; do
    convert_file "$file"
done

echo "Conversion complete!"