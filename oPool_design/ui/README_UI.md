# oPool Design Pipeline - Web UI

A modern, user-friendly web interface for the oPool Design Pipeline, enabling researchers to run antibody library construction workflows through an intuitive browser-based interface.

## Features

### 🚀 **Modern Web Interface**
- Responsive design that works on desktop and mobile devices
- Intuitive step-by-step workflow navigation
- Real-time progress tracking and status updates
- File management and preview capabilities

### 🔬 **Pipeline Steps**
1. **Data Filtering & Extraction** - Upload and filter antibody sequences
2. **Sequence Segments Iteration** - Generate random codon variants
3. **CD-HIT Clustering** - Group similar sequences
4. **Result Selection & Reassembly** - Create assembly libraries
5. **Overlap Region Selection** - Find unique overlap regions
6. **Library Truncation & Primers** - Generate final libraries

### 📁 **File Management**
- Drag-and-drop file upload
- Support for Excel (.xlsx), CSV (.csv), FASTA (.fa), and TSV (.tsv) files
- File preview functionality
- Download capabilities for results
- Organized file storage in uploads/ and results/ folders

### ⚙️ **Configuration**
- **Automatic germline database detection** - No more hardcoded paths!
- Customizable filtering parameters
- Adjustable clustering parameters
- Configurable group sizes and counts
- Environment variable support

## Installation

### Prerequisites
- Python 3.9+
- Conda environment with oPool dependencies (see environment.yml)

### Setup
1. **Activate the oPool environment:**
   ```bash
   conda activate oPool
   ```

2. **Install Flask dependencies:**
   ```bash
   cd ui
   pip install -r requirements.txt
   ```

3. **Create necessary directories:**
   ```bash
   # From the oPool_design directory
   mkdir -p uploads results logs
   ```

## Usage

### Starting the Web UI

#### Option 1: From main directory (Recommended)
```bash
# From oPool_design directory
python launch_ui.py
# or
./launch_ui.sh
```

#### Option 2: From UI directory
```bash
cd ui
python start_ui.py
# or
./start.sh
```

The web interface will be available at: `http://localhost:5001`

### Germline Database Configuration

The UI automatically detects your PyIR germline database location. It searches in this order:

1. **Environment Variable** (highest priority):
   ```bash
   export GERMLINE_DB_PATH="/path/to/your/germline/database"
   ```

2. **Auto-detection** (searches common locations):
   - `~/miniconda3/envs/oPool/lib/python3.9/site-packages/crowelab_pyir/data/germlines/Ig/human`
   - `~/miniconda3/envs/Abs/lib/python3.9/site-packages/crowelab_pyir/data/germlines/Ig/human`
   - `~/anaconda3/envs/oPool/lib/python3.9/site-packages/crowelab_pyir/data/germlines/Ig/human`
   - `~/anaconda3/envs/Abs/lib/python3.9/site-packages/crowelab_pyir/data/germlines/Ig/human`
   - System-wide paths: `/usr/local/share/pyir/germlines/Ig/human`
   - Local project path: `./germlines/Ig/human`

3. **Manual Input**: You can still manually specify the path in the web interface

### Workflow

#### Step 1: Data Upload & Filtering
1. Upload your antibody sequence file (Excel/CSV format)
2. Configure V and D gene family filters
3. The germline database path will be auto-detected
4. Run extraction to get filtered sequences

#### Step 2: Sequence Iteration
1. Select input file from previous step
2. Set pool size (number of random sequences)
3. Choose negative control file
4. Generate codon variants and segments

#### Step 3: CD-HIT Clustering
1. Run CD-HIT clustering on generated segments
2. Monitor progress through log files
3. Results stored in cdhit/ directory

#### Step 4: Result Selection
1. Configure group parameters
2. Select sequences for assembly libraries
3. Generate reassembled antibodies

#### Step 5: Overlap Selection
1. Generate overlap primers around CDR regions
2. Use BLAST+ to find unique overlap regions
3. Select optimal truncation sites

#### Step 6: Final Library
1. Generate final antibody library
2. Add replication primers
3. Ready for DNA synthesis

## File Structure

```
oPool_design/
├── ui/                          # Web UI files
│   ├── app.py                   # Flask web application
│   ├── config.py                # Configuration settings
│   ├── requirements.txt         # Python dependencies
│   ├── templates/               # HTML templates
│   │   └── index.html          # Main UI template
│   ├── start_ui.py             # Python startup script
│   ├── start.sh                # Shell startup script
│   └── README_UI.md            # This file
├── launch_ui.py                 # Main launcher script
├── launch_ui.sh                 # Main shell launcher
├── uploads/                     # User uploaded files
├── results/                     # Pipeline output files
├── script/                      # Original pipeline scripts
├── data/                        # Input data files
└── result/                      # Original pipeline results
```

## Configuration

### Environment Variables
- `FLASK_ENV`: Set to 'development', 'production', or 'testing'
- `FLASK_PORT`: Custom port (default: 5001)
- `FLASK_HOST`: Custom host (default: 0.0.0.0)
- `GERMLINE_DB_PATH`: Custom germline database path
- `SECRET_KEY`: Custom secret key for production

### Default Settings
- **Port**: 5001 (avoiding macOS AirPlay conflict)
- **Upload Limit**: 100MB
- **Pool Size**: 2,000,000 sequences
- **Group Size**: 25 sequences per group
- **Number of Groups**: 12
- **Negative Controls**: 2 per group

## Troubleshooting

### Common Issues

1. **File Upload Fails**
   - Check file format (supported: .xlsx, .csv, .fa, .fasta, .tsv)
   - Ensure file size is reasonable
   - Check uploads/ directory permissions

2. **Pipeline Steps Fail**
   - Verify all dependencies are installed
   - Check script/ directory contains original pipeline scripts
   - Review error messages in the web interface
   - Check console output for detailed errors

3. **Germline Database Not Found**
   - Set `GERMLINE_DB_PATH` environment variable
   - Install PyIR: `conda install -c bioconda pyir`
   - Check if the database exists in common locations
   - Use the auto-detect button in the web interface

4. **CD-HIT Takes Too Long**
   - Reduce dataset size for testing
   - Check cd-hit.sh script parameters
   - Monitor system resources

### Logs
- Pipeline execution logs are displayed in the web interface
- Check the browser console for JavaScript errors
- Flask application logs appear in the terminal

## Development

### Adding New Pipeline Steps
1. Add new route in `ui/app.py`
2. Update HTML template in `ui/templates/index.html`
3. Add corresponding JavaScript functions
4. Test with sample data

### Customizing the UI
- Modify CSS in the `<style>` section of `index.html`
- Update JavaScript functions for new functionality
- Add new Bootstrap components as needed

### Testing
```bash
# Create demo data
cd ui
python demo.py

# Start the server
python app.py
```

## Security Notes

- The web interface runs on localhost by default
- File uploads are restricted to specific extensions
- Input validation is performed on both client and server side
- Consider adding authentication for production use

## Support

For issues with the web interface:
1. Check the troubleshooting section above
2. Review Flask application logs
3. Verify all dependencies are correctly installed
4. Ensure original pipeline scripts are functional

For issues with the underlying pipeline:
- Refer to the main README.md
- Check script documentation
- Verify environment.yml dependencies

## Portability

This UI is designed to work on any system that clones the repository:
- ✅ **No hardcoded paths** - automatically detects your system
- ✅ **Environment variable support** - customize for your setup
- ✅ **Cross-platform compatibility** - works on macOS, Linux, Windows
- ✅ **Auto-detection** - finds PyIR database automatically
- ✅ **Flexible configuration** - easy to adapt to different environments 