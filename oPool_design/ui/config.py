"""
Configuration file for the oPool Design Pipeline Web UI
"""

import os
from pathlib import Path

# Base directory
BASE_DIR = Path(__file__).parent.parent

# Flask configuration
class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'your-secret-key-change-in-production'
    DEBUG = os.environ.get('FLASK_DEBUG', 'True').lower() == 'true'
    PORT = int(os.environ.get('FLASK_PORT', 5001))
    HOST = os.environ.get('FLASK_HOST', '0.0.0.0')
    
    # File upload configuration
    MAX_CONTENT_LENGTH = 100 * 1024 * 1024  # 100MB max file size
    UPLOAD_FOLDER = BASE_DIR / 'uploads'
    RESULT_FOLDER = BASE_DIR / 'ui_results'
    
    # Allowed file extensions
    ALLOWED_EXTENSIONS = {
        'xlsx', 'csv', 'fa', 'fasta', 'tsv', 'txt'
    }
    
    # Pipeline configuration
    DEFAULT_POOL_SIZE = 2000000
    DEFAULT_GROUP_SIZE = 25
    DEFAULT_NUM_GROUPS = 12
    DEFAULT_NUM_NEGATIVE = 2
    
    # CD-HIT configuration
    CDHIT_SIMILARITY_THRESHOLD = 0.85
    CDHIT_WORD_LENGTH = 5
    
    # BLAST configuration
    BLAST_EVALUE = '1e-1'
    BLAST_NUM_THREADS = 8
    BLAST_MAX_HSPS = 2
    
    # Overlap configuration
    OVERLAP_WINDOW_SIZE = 30
    OVERLAP_SLIDING_STEP = 5
    
    # Germline database paths - portable and configurable
    # Users can set GERMLINE_DB_PATH environment variable or use auto-detection
    COMMON_GERMLINE_PATHS = [
        # Environment variable (highest priority)
        os.environ.get('GERMLINE_DB_PATH', ''),
        # Common conda environment paths
        str(Path.home() / 'miniconda3' / 'envs' / 'oPool' / 'lib' / 'python3.9' / 'site-packages' / 'crowelab_pyir' / 'data' / 'germlines' / 'Ig' / 'human'),
        str(Path.home() / 'miniconda3' / 'envs' / 'Abs' / 'lib' / 'python3.9' / 'site-packages' / 'crowelab_pyir' / 'data' / 'germlines' / 'Ig' / 'human'),
        str(Path.home() / 'anaconda3' / 'envs' / 'oPool' / 'lib' / 'python3.9' / 'site-packages' / 'crowelab_pyir' / 'data' / 'germlines' / 'Ig' / 'human'),
        str(Path.home() / 'anaconda3' / 'envs' / 'Abs' / 'lib' / 'python3.9' / 'site-packages' / 'crowelab_pyir' / 'data' / 'germlines' / 'Ig' / 'human'),
        # Dynamic conda environment path detection
        str(Path(os.environ.get('CONDA_PREFIX', '')) / 'lib' / 'python3.9' / 'site-packages' / 'crowelab_pyir' / 'data' / 'germlines' / 'Ig' / 'human'),
        str(Path(os.environ.get('CONDA_PREFIX', '')) / 'lib' / 'python3.8' / 'site-packages' / 'crowelab_pyir' / 'data' / 'germlines' / 'Ig' / 'human'),
        str(Path(os.environ.get('CONDA_PREFIX', '')) / 'lib' / 'python3.7' / 'site-packages' / 'crowelab_pyir' / 'data' / 'germlines' / 'Ig' / 'human'),
        # System-wide paths
        '/usr/local/share/pyir/germlines/Ig/human',
        '/opt/pyir/germlines/Ig/human',
        # Local project path (for development)
        str(BASE_DIR / 'germlines' / 'Ig' / 'human'),
        # Empty string as fallback
        ''
    ]
    
    # Default V and D gene families to remove
    DEFAULT_V_GENE_FAMILIES = ['IGHV1-69', 'IGHV6-1', 'IGHV1-18']
    DEFAULT_D_GENE_FAMILIES = ['IGHD3-9']
    
    # Clonotype filtering configuration
    DEFAULT_SKIP_CLONOTYPE_FILTER = False  # Set to True to skip clonotype filtering by default
    
    # Logging configuration
    LOG_LEVEL = os.environ.get('LOG_LEVEL', 'INFO')
    LOG_FILE = BASE_DIR / 'logs' / 'opool_ui.log'
    
    # Session configuration
    PERMANENT_SESSION_LIFETIME = 3600  # 1 hour
    
    # Security configuration
    SESSION_COOKIE_SECURE = False  # Set to True in production with HTTPS
    SESSION_COOKIE_HTTPONLY = True
    SESSION_COOKIE_SAMESITE = 'Lax'

class DevelopmentConfig(Config):
    DEBUG = True
    TESTING = False

class ProductionConfig(Config):
    DEBUG = False
    TESTING = False
    SECRET_KEY = os.environ.get('SECRET_KEY')
    SESSION_COOKIE_SECURE = True

class TestingConfig(Config):
    TESTING = True
    DEBUG = True
    WTF_CSRF_ENABLED = False

# Configuration dictionary
config = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
    'testing': TestingConfig,
    'default': DevelopmentConfig
}

def get_config():
    """Get configuration based on environment"""
    config_name = os.environ.get('FLASK_ENV', 'default')
    return config.get(config_name, config['default'])

def get_germline_path():
    """Get the first available germline database path"""
    config = get_config()
    
    # Try to find PyIR installation and derive germline path
    try:
        import pyir
        pyir_path = Path(pyir.__file__).parent
        # Look for germlines in the PyIR package directory
        germline_paths = [
            pyir_path / 'data' / 'germlines' / 'Ig' / 'human',
            pyir_path.parent / 'data' / 'germlines' / 'Ig' / 'human',
        ]
        
        for path in germline_paths:
            if path.exists():
                return str(path)
    except ImportError:
        pass
    
    # Filter out empty paths and check if directories exist
    valid_paths = []
    for path in config.COMMON_GERMLINE_PATHS:
        if path and os.path.exists(path):
            valid_paths.append(path)
    
    if valid_paths:
        return valid_paths[0]
    
    # If no valid paths found, return None
    return None

def create_directories():
    """Create necessary directories"""
    directories = [
        Config.UPLOAD_FOLDER,
        Config.RESULT_FOLDER,
        Config.LOG_FILE.parent
    ]
    
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)

def validate_config():
    """Validate configuration settings"""
    errors = []
    
    # Check if upload and result folders are writable
    try:
        Config.UPLOAD_FOLDER.mkdir(parents=True, exist_ok=True)
        test_file = Config.UPLOAD_FOLDER / 'test.txt'
        test_file.write_text('test')
        test_file.unlink()
    except Exception as e:
        errors.append(f"Upload folder not writable: {e}")
    
    # Check if germline database is accessible
    germline_path = get_germline_path()
    if not germline_path:
        errors.append("Germline database not found. Please set GERMLINE_DB_PATH environment variable or ensure PyIR is properly installed.")
    
    return errors

def get_environment_setup_instructions():
    """Get instructions for setting up the environment"""
    instructions = []
    
    # Check if we're in a conda environment
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        instructions.append(f"✅ Conda environment detected: {conda_prefix}")
        
        # Try to find PyIR in the current environment
        try:
            import pyir
            pyir_path = Path(pyir.__file__).parent
            germline_path = pyir_path / 'data' / 'germlines' / 'Ig' / 'human'
            
            if germline_path.exists():
                instructions.append(f"✅ PyIR germline database found at: {germline_path}")
                instructions.append(f"   You can set: export GERMLINE_DB_PATH='{germline_path}'")
            else:
                instructions.append(f"⚠️  PyIR installed but germline database not found at: {germline_path}")
        except ImportError:
            instructions.append("❌ PyIR not installed in current environment")
            instructions.append("   Install with: pip install crowelab_pyir")
    else:
        instructions.append("⚠️  Not in a conda environment")
        instructions.append("   Consider activating your conda environment first")
    
    return instructions
    
    try:
        Config.RESULT_FOLDER.mkdir(parents=True, exist_ok=True)
        test_file = Config.RESULT_FOLDER / 'test.txt'
        test_file.write_text('test')
        test_file.unlink()
    except Exception as e:
        errors.append(f"Result folder not writable: {e}")
    
    # Check if pipeline scripts exist
    script_dir = BASE_DIR / 'script'
    required_scripts = [
        'extract.py',
        'iteration.py',
        'cd-hit.sh',
        'cdhit_result.py',
        'Overlap_check.py',
        'ChunkByOverlap.py'
    ]
    
    for script in required_scripts:
        if not (script_dir / script).exists():
            errors.append(f"Missing pipeline script: {script}")
    
    # Check germline database
    germline_path = get_germline_path()
    if not germline_path:
        errors.append("No valid germline database path found. Please set GERMLINE_DB_PATH environment variable or install PyIR.")
    
    return errors 