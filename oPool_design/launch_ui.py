#!/usr/bin/env python3
"""
Launcher script for the oPool Design Pipeline Web UI
This script can be run from the oPool_design directory
"""

import os
import sys
import subprocess
from pathlib import Path

def main():
    """Main launcher function"""
    print("🔬 oPool Design Pipeline Web UI Launcher")
    print("=" * 45)
    
    # Get the directory where this script is located
    script_dir = Path(__file__).parent
    ui_dir = script_dir / 'ui'
    
    # Check if we're in the right directory
    if not ui_dir.exists():
        print("❌ Error: ui directory not found!")
        print(f"Expected path: {ui_dir}")
        sys.exit(1)
    
    if not (ui_dir / 'app.py').exists():
        print("❌ Error: app.py not found in ui directory!")
        print(f"Expected path: {ui_dir / 'app.py'}")
        sys.exit(1)
    
    print(f"✅ Found UI directory: {ui_dir}")
    print(f"📁 Current working directory: {os.getcwd()}")
    
    # Change to ui directory and start the application
    print("\n🚀 Launching Web UI...")
    os.chdir(ui_dir)
    
    try:
        # Import and run the Flask app
        from app import app
        print("✅ Flask app imported successfully")
        print("🌐 Starting server on http://localhost:5001")
        print("📱 The interface will open in your default browser")
        print("⏹️  Press Ctrl+C to stop the server\n")
        
        app.run(debug=True, host='0.0.0.0', port=5001)
        
    except ImportError as e:
        print(f"❌ Error importing Flask app: {e}")
        print("Please install dependencies: pip install -r ui/requirements.txt")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error starting server: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main() 