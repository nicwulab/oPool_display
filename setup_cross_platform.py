#!/usr/bin/env python3
"""
Cross-platform setup script for oPool Design Pipeline
Supports both Linux and macOS
"""

import os
import sys
import subprocess
import platform
from pathlib import Path

def print_header():
    """Print setup header"""
    print("🚀 oPool Design Pipeline - Cross-Platform Setup")
    print("=" * 50)
    print(f"🖥️  Platform: {platform.system()} {platform.release()}")
    print(f"🐍 Python: {sys.version}")
    print()

def check_conda():
    """Check if conda is available"""
    try:
        result = subprocess.run(['conda', '--version'], 
                              capture_output=True, text=True, check=True)
        print(f"✅ Conda found: {result.stdout.strip()}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("❌ Conda not found. Please install Miniconda or Anaconda first.")
        print("💡 Download from: https://docs.conda.io/en/latest/miniconda.html")
        return False

def create_conda_env():
    """Create conda environment"""
    env_name = "oPool"
    
    # Check if environment already exists
    try:
        result = subprocess.run(['conda', 'env', 'list'], 
                              capture_output=True, text=True, check=True)
        if env_name in result.stdout:
            print(f"⚠️  Environment '{env_name}' already exists")
            response = input("Do you want to recreate it? (y/N): ").strip().lower()
            if response == 'y':
                print(f"🗑️  Removing existing environment '{env_name}'...")
                subprocess.run(['conda', 'env', 'remove', '-n', env_name, '-y'], check=True)
            else:
                print(f"✅ Using existing environment '{env_name}'")
                return True
    except subprocess.CalledProcessError:
        print("⚠️  Could not check existing environments")
    
    print(f"🔧 Creating conda environment '{env_name}'...")
    
    try:
        # Create environment with Python 3.9
        subprocess.run(['conda', 'create', '-n', env_name, 'python=3.9', '-y'], check=True)
        print(f"✅ Environment '{env_name}' created successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to create environment: {e}")
        return False

def install_dependencies():
    """Install dependencies from environment.yml"""
    print("📦 Installing dependencies...")
    
    try:
        # Install from environment.yml
        subprocess.run(['conda', 'env', 'update', '-f', 'environment.yml'], check=True)
        print("✅ Dependencies installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install dependencies: {e}")
        return False

def get_germline_path():
    """Get the path to the PyIR germline database"""
    try:
        import crowelab_pyir
        pyir_path = os.path.dirname(crowelab_pyir.__file__)
        germline_path = os.path.join(pyir_path, 'data', 'germlines')
        
        if os.path.exists(germline_path):
            return germline_path
        else:
            return None
    except ImportError:
        return None

def check_germline_completeness(germline_path):
    """Check if germline database is complete enough to use"""
    if not germline_path or not os.path.exists(germline_path):
        return False

    # Check for common species directories
    species_dirs = ['human', 'mouse', 'rat']
    available_species = []
    
    for species in species_dirs:
        species_path = os.path.join(germline_path, 'Ig', species)
        if os.path.exists(species_path):
            available_species.append(species)
    
    if available_species:
        print(f"🔍 Found germline data for: {', '.join(available_species)}")
    else:
        print("🔍 Found germline data for: none")
    
    # Consider it usable if we have at least human data
    return 'human' in available_species

def setup_pyir():
    """Setup PyIR with better error handling"""
    print("🔧 Setting up PyIR...")
    
    try:
        # Install crowelab-pyir via pip (it's not in conda)
        subprocess.run(['pip', 'install', 'crowelab-pyir'], check=True)
        print("✅ PyIR installed successfully")
        
        # Setup germline database with timeout handling
        print("🔍 Setting up germline database...")
        print("⚠️  Note: IMGT.org downloads can be slow - this may take 20+ minutes")
        
        try:
            # Try with longer timeout
            result = subprocess.run(['pyir', 'setup'], 
                                  capture_output=True, text=True, 
                                  timeout=3600)  # 1 hour timeout
            
            if result.returncode == 0:
                print("✅ Germline database setup completed successfully")
                return True
            else:
                print(f"⚠️  PyIR setup returned code {result.returncode}")
                print("⚠️  Setup may be incomplete - checking what's available...")
                
                # Check if we can proceed with partial data
                germline_path = get_germline_path()
                if germline_path and check_germline_completeness(germline_path):
                    print("✅ Partial setup detected - proceeding with available data")
                    return True
                else:
                    print("❌ Setup failed and no usable database found")
                    return False
                    
        except subprocess.TimeoutExpired:
            print("⏰ PyIR setup timed out (1 hour)")
            print("🔄 Checking for partial setup...")
            
            # Check if we can proceed with what was downloaded
            germline_path = get_germline_path()
            if germline_path and check_germline_completeness(germline_path):
                print("✅ Partial setup detected - proceeding with available data")
                return True
            else:
                print("❌ Setup failed and database is not usable")
                print("💡 You can try running 'pyir setup' manually later")
                return False
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to setup PyIR: {e}")
        return False

def create_directories():
    """Create necessary directories"""
    print("🔍 Creating directories...")
    
    dirs = [
        'oPool_design/uploads',
        'oPool_design/ui_results', 
        'oPool_design/logs'
    ]
    
    for dir_path in dirs:
        if os.path.exists(dir_path):
            print(f"🔍 Exists: {dir_path}")
        else:
            os.makedirs(dir_path, exist_ok=True)
            print(f"✅ Created: {dir_path}")

def check_dependencies():
    """Check if all required packages are installed"""
    print("🔍 Checking dependencies...")
    
    required_packages = [
        'flask', 'pandas', 'numpy', 'matplotlib', 'seaborn',
        'biopython', 'cd-hit', 'blast', 'hmmer'
    ]
    
    missing = []
    for package in required_packages:
        try:
            __import__(package.replace('-', '_'))
            print(f"✅ {package}")
        except ImportError:
            missing.append(package)
            print(f"❌ {package}")
    
    if missing:
        missing_str = ', '.join(missing)
        print(f"\n❌ Missing packages: {missing_str}")
        print("�� Install missing packages with: pip install " + " ".join(missing))
        return False
    else:
        print("✅ All required packages are available")
        return True

def main():
    """Main setup function"""
    print_header()
    
    # Check conda
    if not check_conda():
        return False
    
    # Create environment
    if not create_conda_env():
        return False
    
    # Install dependencies
    if not install_dependencies():
        return False
    
    # Setup PyIR
    if not setup_pyir():
        print("⚠️  PyIR setup failed - you may need to run 'pyir setup' manually")
    
    # Create directories
    create_directories()
    
    # Check dependencies
    check_dependencies()
    
    print("\n🎉 Setup completed!")
    print("\n📋 Next steps:")
    print("1. Activate environment: conda activate oPool")
    print("2. Navigate to oPool_design: cd oPool_design")
    print("3. Launch UI: python launch_ui.py")
    print("\n💡 If you encounter issues, check the troubleshooting guide in CROSS_PLATFORM_SETUP.md")
    
    return True

if __name__ == "__main__":
    try:
        success = main()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n\n⚠️  Setup interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        sys.exit(1)
