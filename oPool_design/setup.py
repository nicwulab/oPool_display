#!/usr/bin/env python3
"""
Setup script for oPool Design Pipeline
This script redirects to the parent directory setup for cross-platform compatibility
"""

import os
import sys
from pathlib import Path

def main():
    """Redirect to parent directory setup"""
    print("🚀 oPool Design Pipeline Setup")
    print("=" * 40)
    print()
    print("📁 This setup script has been moved to the parent directory for better organization.")
    print()
    print("🔧 To set up the environment, please run from the parent directory:")
    print()
    print("   cd ..")
    print("   python setup_cross_platform.py")
    print("   # OR")
    print("   ./setup_cross_platform.sh")
    print()
    print("📋 After environment setup, return here to launch the UI:")
    print("   cd oPool_design")
    print("   python launch_ui.py")
    print()
    print("📖 For detailed setup instructions, see:")
    print("   ../CROSS_PLATFORM_SETUP.md")
    print()
    print("🔗 Or visit the main README:")
    print("   ../README.md")

if __name__ == "__main__":
    main() 