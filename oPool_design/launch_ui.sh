#!/bin/bash

# oPool Design Pipeline Web UI Launcher
# This script can be run from the oPool_design directory

echo "🔬 oPool Design Pipeline Web UI Launcher"
echo "========================================"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UI_DIR="$SCRIPT_DIR/ui"

# Check if we're in the right directory
if [ ! -d "$UI_DIR" ]; then
    echo "❌ Error: ui directory not found!"
    echo "Expected path: $UI_DIR"
    exit 1
fi

if [ ! -f "$UI_DIR/app.py" ]; then
    echo "❌ Error: app.py not found in ui directory!"
    echo "Expected path: $UI_DIR/app.py"
    exit 1
fi

echo "✅ Found UI directory: $UI_DIR"
echo "📁 Current working directory: $(pwd)"

# Change to ui directory and start the application
echo ""
echo "🚀 Launching Web UI..."
cd "$UI_DIR"

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo "❌ Error: Python not found!"
    echo "Please install Python 3.9+ or activate your conda environment"
    exit 1
fi

# Check if conda environment is activated
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "⚠️  Warning: No conda environment detected"
    echo "Consider activating the oPool environment: conda activate oPool"
else
    echo "✅ Conda environment: $CONDA_DEFAULT_ENV"
fi

# Check dependencies
echo "🔍 Checking dependencies..."
python -c "import flask, pandas, numpy, openpyxl" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "❌ Missing dependencies. Installing..."
    pip install -r requirements.txt
fi

# Start the application
echo "🚀 Starting web UI..."
echo "📱 The interface will open in your default browser"
echo "🌐 Server will be available at: http://localhost:5001"
echo "⏹️  Press Ctrl+C to stop the server"
echo ""

python app.py 