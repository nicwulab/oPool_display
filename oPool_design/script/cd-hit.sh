#!/bin/bash

# CD-HIT script for organized pipeline structure
# This script can be run from any directory

# Parse command line arguments
MIN_THRESHOLD=${1:-0.6}
MAX_THRESHOLD=${2:-0.85}
INCREMENT=${3:-0.05}

echo "CD-HIT Parameters:"
echo "  Min Threshold: $MIN_THRESHOLD"
echo "  Max Threshold: $MAX_THRESHOLD"
echo "  Increment: $INCREMENT"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Find the most recent iteration output file from step2
INPUT_FILE=$(find "$PROJECT_ROOT/ui_results/step2/" -name "*.fa" -type f | head -1)

if [ -z "$INPUT_FILE" ]; then
    echo "Error: No FASTA input file found in ui_results/step2/"
    echo "Please complete Step 2 (Iteration) first."
    exit 1
fi

echo "Using input file: $INPUT_FILE"

# Create step3 directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/ui_results/step3"

# Change to step3 directory for output
cd "$PROJECT_ROOT/ui_results/step3"

# Generate threshold values using awk for floating point arithmetic
echo "Generating similarity thresholds from $MIN_THRESHOLD to $MAX_THRESHOLD with increment $INCREMENT"

# Use awk to generate the sequence of thresholds
awk -v min="$MIN_THRESHOLD" -v max="$MAX_THRESHOLD" -v inc="$INCREMENT" '
BEGIN {
    for (x = min; x <= max + inc/2; x += inc) {
        # Format to remove trailing zeros for compatibility with cdhit_result_modified.py
        if (x == int(x)) {
            printf "%.0f\n", x
        } else {
            val = sprintf("%.2f", x)
            gsub(/\.?0+$/, "", val)  # Remove trailing zeros
            print val
        }
    }
}' | while read x; do
    echo "Running CD-HIT with similarity threshold: $x"
    cd-hit -i "$INPUT_FILE" -o "HA_Abs_paired_$x.fa" -c $x -M 50000 -d 0 -T 60 -n 3.5 -aL 0.9 -s 0.95 -uS 0.2 -sc 1 -sf 1
    
    if [ $? -ne 0 ]; then
        echo "Error: CD-HIT failed for threshold $x"
        exit 1
    fi
done

echo "CD-HIT clustering completed. Results saved in ui_results/step3/"
