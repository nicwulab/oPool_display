#!/usr/bin/env python3
"""
Setup script for the oPool Design Pipeline Web UI
This script helps configure the UI for your system
"""

import os
import sys
import subprocess
from pathlib import Path

def check_pyir_installation():
    """Check if PyIR is installed and find the germline database"""
    print("🔍 Checking PyIR installation...")
    
    try:
        # Try to import pyir
        import pyir
        print("✅ PyIR is installed")
        
        # Try to find the germline database
        from config import get_germline_path
        germline_path = get_germline_path()
        
        if germline_path:
            print(f"✅ Germline database found at: {germline_path}")
            return True
        else:
            print("⚠️  Germline database not found automatically")
            return False
            
    except ImportError:
        print("❌ PyIR is not installed")
        print("Please install it using: conda install -c bioconda pyir")
        return False

def check_dependencies():
    """Check if required dependencies are installed"""
    print("\n🔍 Checking Python dependencies...")
    
    required_packages = ['flask', 'pandas', 'numpy', 'openpyxl']
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"✅ {package}")
        except ImportError:
            print(f"❌ {package}")
            missing_packages.append(package)
    
    if missing_packages:
        print(f"\n❌ Missing packages: {', '.join(missing_packages)}")
        print("Installing missing packages...")
        
        try:
            subprocess.run([sys.executable, '-m', 'pip', 'install'] + missing_packages, check=True)
            print("✅ Dependencies installed successfully")
        except subprocess.CalledProcessError:
            print("❌ Failed to install dependencies")
            return False
    
    return True

def create_directories():
    """Create necessary directories"""
    print("\n📁 Creating directories...")
    
    directories = ['uploads', 'results', 'logs']
    
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"✅ Created: {directory}")
        else:
            print(f"📁 Exists: {directory}")

def check_pipeline_scripts():
    """Check if pipeline scripts exist"""
    print("\n🔍 Checking pipeline scripts...")
    
    script_dir = Path('script')
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
        script_path = script_dir / script
        if not script_path.exists():
            missing_scripts.append(script)
            print(f"❌ Missing: {script}")
        else:
            print(f"✅ Found: {script}")
    
    if missing_scripts:
        print(f"\n⚠️  Missing scripts: {', '.join(missing_scripts)}")
        print("Make sure you're running this from the oPool_design directory")
        return False
    
    return True

def setup_environment():
    """Set up environment variables"""
    print("\n⚙️  Setting up environment...")
    
    # Check if conda environment is activated
    conda_env = os.environ.get('CONDA_DEFAULT_ENV')
    if conda_env:
        print(f"✅ Conda environment: {conda_env}")
    else:
        print("⚠️  No conda environment detected")
        print("Consider activating the oPool environment: conda activate oPool")
    
    # Check for custom germline path
    germline_path = os.environ.get('GERMLINE_DB_PATH')
    if germline_path:
        print(f"✅ Custom germline path: {germline_path}")
    else:
        print("ℹ️  No custom germline path set (will use auto-detection)")
    
    # Suggest environment variable setup
    print("\n💡 To set a custom germline database path, add this to your shell profile:")
    print("export GERMLINE_DB_PATH=\"/path/to/your/germline/database\"")

def main():
    """Main setup function"""
    print("🔬 oPool Design Pipeline Web UI Setup")
    print("=" * 45)
    
    # Check if we're in the right directory
    if not os.path.exists('ui/app.py'):
        print("❌ Error: ui/app.py not found!")
        print("Please run this script from the oPool_design directory")
        sys.exit(1)
    
    # Run all checks
    success = True
    
    if not check_dependencies():
        success = False
    
    if not check_pipeline_scripts():
        success = False
    
    check_pyir_installation()
    
    create_directories()
    setup_environment()
    
    if success:
        print("\n🎉 Setup completed successfully!")
        print("\n🚀 To start the web UI, run:")
        print("  python launch_ui.py")
        print("  or")
        print("  ./launch_ui.sh")
        print("\n🌐 The interface will be available at: http://localhost:5001")
    else:
        print("\n⚠️  Setup completed with warnings")
        print("Please resolve the issues above before starting the UI")
        sys.exit(1)

if __name__ == '__main__':
    main() 