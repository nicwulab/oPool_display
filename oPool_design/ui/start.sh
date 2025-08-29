#!/bin/bash

# oPool Design Pipeline Web UI Startup Script

echo "🔬 oPool Design Pipeline Web UI"
echo "================================"

# Check if we're in the right directory
if [ ! -f "app.py" ]; then
    echo "❌ Error: app.py not found!"
    echo "Please run this script from the ui directory"
    exit 1
fi

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

# Create necessary directories in parent directory
echo "📁 Creating directories..."
PARENT_DIR="$(dirname "$(pwd)")"
mkdir -p "$PARENT_DIR/uploads" "$PARENT_DIR/results" "$PARENT_DIR/logs"

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