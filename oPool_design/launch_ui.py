#!/usr/bin/env python3
"""
oPool Design Pipeline Web UI Launcher
Launches the Flask web interface for the pipeline
"""

import os
import sys
from pathlib import Path

def main():
    """Main launcher function"""
    print("🔬 oPool Design Pipeline Web UI Launcher")
    print("=" * 50)
    
    # Get the current working directory
    current_dir = Path.cwd()
    print(f"📁 Current working directory: {current_dir}")
    
    # Check if we're in the right directory
    ui_dir = current_dir / "ui"
    if not ui_dir.exists():
        print("❌ UI directory not found. Please run this script from the oPool_design directory.")
        return False
    
    print(f"✅ Found UI directory: {ui_dir}")
    
    # Check if app.py exists
    app_file = ui_dir / "app.py"
    if not app_file.exists():
        print("❌ app.py not found in UI directory.")
        return False
    
    print("🚀 Launching Web UI...")
    
    try:
        # Add ui directory to Python path
        sys.path.insert(0, str(ui_dir))
        
        # Import and run the Flask app
        from app import app
        
        print("✅ Flask app imported successfully")
        print("🌐 Starting web server...")
        print("🔗 Open your browser to: http://127.0.0.1:5001")
        print("⏹️  Press Ctrl+C to stop the server")
        
        # Run the Flask app
        app.run(
            debug=False,
            host='127.0.0.1',
            port=5001,
            use_reloader=False
        )
        
    except ImportError as e:
        print(f"❌ Error importing Flask app: {e}")
        print("Please check that all dependencies are installed")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False

if __name__ == "__main__":
    try:
        success = main()
        if not success:
            sys.exit(1)
    except KeyboardInterrupt:
        print("\n\n⚠️  Server stopped by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        sys.exit(1) 