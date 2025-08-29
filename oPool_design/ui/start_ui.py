#!/usr/bin/env python3
"""
Startup script for the oPool Design Pipeline Web UI
"""

import os
import sys
import subprocess
import webbrowser
import time

def check_dependencies():
    """Check if required dependencies are installed"""
    required_packages = ['flask', 'pandas', 'numpy', 'openpyxl']
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print(f"❌ Missing required packages: {', '.join(missing_packages)}")
        print("Please install them using: pip install -r requirements.txt")
        return False
    
    print("✅ All required packages are installed")
    return True

def check_directories():
    """Create necessary directories"""
    # Get the parent directory (oPool_design)
    parent_dir = os.path.dirname(os.path.abspath(__file__))
    directories = [
        os.path.join(parent_dir, 'uploads'),
        os.path.join(parent_dir, 'results'),
        os.path.join(parent_dir, 'logs')
    ]
    
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"📁 Created directory: {directory}")
        else:
            print(f"📁 Directory exists: {directory}")

def check_scripts():
    """Check if pipeline scripts exist"""
    # Get the parent directory (oPool_design)
    parent_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = os.path.join(parent_dir, 'script')
    required_scripts = [
        'extract.py',
        'iteration.py',
        'cd-hit.sh',
        'cdhit_result.py',
        'Overlap_check.py',
        'ChunkByOverlap.py'
    ]
    
    missing_scripts = []
    for script in required_scripts:
        script_path = os.path.join(script_dir, script)
        if not os.path.exists(script_path):
            missing_scripts.append(script)
    
    if missing_scripts:
        print(f"⚠️  Missing pipeline scripts: {', '.join(missing_scripts)}")
        print("Make sure you're running this from the oPool_design directory")
        return False
    
    print("✅ All pipeline scripts are present")
    return True

def start_server():
    """Start the Flask web server"""
    print("\n🚀 Starting oPool Design Pipeline Web UI...")
    print("📱 The interface will open in your default browser")
    print("🌐 Server will be available at: http://localhost:5001")
    print("⏹️  Press Ctrl+C to stop the server\n")
    
    # Wait a moment for server to start
    time.sleep(2)
    
    # Open browser
    try:
        webbrowser.open('http://localhost:5001')
    except:
        print("⚠️  Could not open browser automatically")
        print("Please open http://localhost:5001 in your browser")
    
    # Start Flask app
    try:
        from app import app
        app.run(debug=True, host='0.0.0.0', port=5001)
    except KeyboardInterrupt:
        print("\n🛑 Server stopped by user")
    except Exception as e:
        print(f"❌ Error starting server: {e}")

def main():
    """Main startup function"""
    print("🔬 oPool Design Pipeline Web UI")
    print("=" * 40)
    
    # Check if we're in the right directory
    if not os.path.exists('app.py'):
        print("❌ Error: app.py not found!")
        print("Please run this script from the ui directory")
        sys.exit(1)
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Check and create directories
    check_directories()
    
    # Check pipeline scripts
    if not check_scripts():
        print("⚠️  Continuing anyway, but some features may not work")
    
    # Start the server
    start_server()

if __name__ == '__main__':
    main() 