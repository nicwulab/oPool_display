# Cross-Platform Setup Guide for oPool Design Pipeline

## 🌍 **Overview**

The oPool Design Pipeline is designed to work seamlessly on both **Linux** and **macOS** operating systems. This guide explains how the cross-platform compatibility is achieved and provides setup instructions for both platforms.

## 🏗️ **Architecture**

### **Single Environment File Approach**
- **`environment.yml`**: Unified conda environment specification
- **No platform-specific files**: Single file works on both systems
- **Automatic dependency resolution**: Conda handles platform-specific package versions

### **Cross-Platform Scripts**
- **`setup_cross_platform.py`**: Python-based setup script
- **`setup_cross_platform.py`**: Python-based setup script
- **Automatic platform detection**: Scripts adapt to the current OS

## 🚀 **Quick Start**

### **For All Platforms**
```bash
# Clone the repository
git clone <repository-url>
cd oPool_display

# Run automated setup
./setup_cross_platform.py
# OR
python setup_cross_platform.py
```

### **Detailed Setup Process**

#### **Step 1: Cross-Platform Setup**
Run the cross-platform setup script. If this is the first time, the script will need to be rerun with "No" selected to allow the installation of other dependencies.

```bash
./setup_cross_platform.py
conda activate oPool
./setup_cross_platform.py
```

#### **Step 2: PyIR Configuration**
**PyIR NOTE**: It is important to fix the line 112 in `setup_germline_library.py` from "http" to "https" to download fasta reference sequences

#### **Step 3: UI Dependencies**
```bash
pip install -r ui/requirements.txt
```

#### **Step 4: Launch the Interface**
```bash
cd oPool_design
python launch_ui.py
# Access at: http://127.0.0.1:5001
```

## 📋 **Platform Support Matrix**

| Feature | macOS | Linux | Notes |
|---------|-------|-------|-------|
| **Core Pipeline** | ✅ | ✅ | Identical functionality |
| **Web UI** | ✅ | ✅ | Flask-based, platform-agnostic |
| **PyIR Integration** | ✅ | ✅ | Automatic path detection |
| **BLAST/CD-HIT** | ✅ | ✅ | Conda packages work on both |
| **File Operations** | ✅ | ✅ | Python handles path differences |
| **Environment Management** | ✅ | ✅ | Conda works identically |

## 🔧 **Technical Implementation**

### **1. Environment Management**
```yaml
# environment.yml - Works on both platforms
name: oPool
channels:
  - bioconda      # Cross-platform bioinformatics packages
  - conda-forge   # Modern, well-maintained packages
  - defaults      # Core conda packages
dependencies:
  - python=3.9   # Available on both platforms
  - cd-hit       # Bioconda provides platform-specific binaries
  - blast        # Bioconda provides platform-specific binaries
  - pip:
    - crowelab-pyir  # Pure Python package, platform-agnostic
```

### **2. Path Handling**
```python
# Automatic path detection in config.py
import os
from pathlib import Path

# Uses environment variables instead of hardcoded paths
CONDA_PREFIX = os.environ.get('CONDA_PREFIX', '')
GERMLINE_DB_PATH = os.environ.get('GERMLINE_DB_PATH', '')

# Platform-agnostic path construction
germline_path = Path(CONDA_PREFIX) / 'lib' / 'python3.9' / 'site-packages' / 'crowelab_pyir' / 'data' / 'germlines' / 'Ig' / 'human'
```

### **3. Platform Detection**
```python
import platform

def detect_platform():
    system = platform.system().lower()
    if system == "darwin":
        return "macos"
    elif system == "linux":
        return "linux"
    else:
        return "unknown"
```

## 🐧 **Linux-Specific Considerations**

### **System Dependencies**
Some Linux distributions may require additional system packages:
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install build-essential libssl-dev libffi-dev

# CentOS/RHEL
sudo yum groupinstall "Development Tools"
sudo yum install openssl-devel libffi-devel

# Arch Linux
sudo pacman -S base-devel openssl libffi
```

### **Conda Installation**
```bash
# Download Miniconda for Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install
bash Miniconda3-latest-Linux-x86_64.sh

# Initialize
source ~/.bashrc
```

## 🍎 **macOS-Specific Considerations**

### **System Dependencies**
macOS typically has fewer system dependency requirements:
```bash
# Install Homebrew (if not already installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Xcode Command Line Tools (if needed)
xcode-select --install
```

### **Conda Installation**
```bash
# Download Miniconda for macOS
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

# Install
bash Miniconda3-latest-MacOSX-x86_64.sh

# Initialize
source ~/.zshrc  # or ~/.bash_profile
```

## 🔍 **Troubleshooting**

### **Common Issues**

#### **1. Conda Not Found**
```bash
# Add conda to PATH
export PATH="$HOME/miniconda3/bin:$PATH"

# Or initialize conda
conda init
```

#### **2. Package Installation Failures**
```bash
# Update conda
conda update conda

# Clear package cache
conda clean --all

# Try alternative channels
conda install -c conda-forge <package-name>
```

#### **3. Permission Issues**
```bash
# Check file permissions
ls -la setup_cross_platform.py

# Make executable
chmod +x setup_cross_platform.py
```

### **Platform-Specific Issues**

#### **Linux**
- **GLIBC version**: Ensure GLIBC >= 2.14
- **Library conflicts**: Use conda environments to avoid system library conflicts
- **Memory limits**: Some operations may require increased memory limits

#### **macOS**
- **Security policies**: May need to allow conda in Security & Privacy settings
- **Path issues**: Ensure conda is in PATH for your shell
- **ARM64 support**: M1/M2 Macs are fully supported

## 📚 **Advanced Configuration**

### **Custom Environment Variables**
```bash
# Set custom paths
export GERMLINE_DB_PATH="/custom/path/to/germlines"
export FLASK_PORT=5002
export FLASK_DEBUG=True
```

### **Multiple Python Versions**
```yaml
# environment.yml supports multiple Python versions
dependencies:
  - python>=3.8,<3.10  # Supports Python 3.8 and 3.9
```

### **Custom Channels**
```yaml
# Add custom conda channels if needed
channels:
  - custom-channel
  - bioconda
  - conda-forge
  - defaults
```

## 🧪 **Testing Cross-Platform Compatibility**

### **Automated Testing**
```bash
# Test on different platforms
docker run -it --rm continuumio/miniconda3:latest bash
# Test Linux compatibility

# On macOS, test with different Python versions
conda create -n test-py38 python=3.8
conda create -n test-py39 python=3.9
```

### **Manual Verification**
```bash
# Check platform detection
python -c "import platform; print(platform.system())"

# Verify PyIR installation
python -c "import pyir; print('PyIR works')"

# Test UI launch
cd oPool_design
python launch_ui.py
```

## 📖 **References**

- [Conda Documentation](https://docs.conda.io/)
- [Bioconda Documentation](https://bioconda.github.io/)
- [Cross-Platform Python Development](https://docs.python.org/3/library/platform.html)
- [Pathlib Documentation](https://docs.python.org/3/library/pathlib.html)


**Note**: This setup has been tested on Ubuntu 20.04+, CentOS 7+, macOS 10.15+, and macOS 12+ (including M1/M2 Macs). 