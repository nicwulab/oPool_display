# Setup Details & Test Results

## 🏗️ **Directory Organization**

The oPool Design Pipeline setup has been reorganized for better cross-platform compatibility and cleaner project structure.

### **Repository Root (`oPool_display/`)**
```
oPool_display/
├── environment.yml                    # Cross-platform conda environment
├── setup_cross_platform.py          # Python setup script
├── setup_cross_platform.py          # Python setup script  
├── CROSS_PLATFORM_SETUP.md          # Detailed setup guide
├── README.md                         # Main project README
└── oPool_design/                    # Pipeline implementation
```

### **Pipeline Directory (`oPool_design/`)**
```
oPool_design/
├── ui/                              # Web UI implementation
├── HA_screen/script/                          # Pipeline scripts (extract.py, iteration.py, etc.)
├── HA_screen/data/                            # Test datasets (TableS1.csv, etc.)
├── ui_results/                      # Pipeline outputs
├── setup.py                         # Setup redirect script
├── launch_ui.py                     # UI launcher
└── README.md                        # Pipeline-specific README
```

## 🚀 **Setup Process**

### **Step 1: Environment Setup (from repository root)**
Please follow environment_setup.md

### **Step 2: Launch UI (from pipeline directory)**
```bash
# Navigate to pipeline directory
cd oPool_design

# Launch the Web UI
python launch_ui.py
```


### **UI Features:**
- **Step-by-step workflow navigation**: Guided process through all 6 pipeline steps
- **File management**: Upload, preview, and download capabilities
- **Real-time progress tracking**: Live updates during processing
- **Parameter configuration**: Customizable settings for each step
- **Optional clonotype filtering**: Toggle for 3x more sequence retention
- **Cross-platform compatibility**: Works on macOS and Linux

## 📋 **Input File Requirements**

### **Supported File Formats**
- **Excel files**: `.xlsx` format (recommended for initial data)
- **CSV files**: `.csv` format with comma separation
- **TSV files**: `.tsv` format with tab separation
- **FASTA files**: `.fa` or `.fasta` format for sequence data

### **Input Data Structure**
The pipeline expects antibody sequence data with the following columns:

| Column | Description | Required | Example |
|--------|-------------|----------|---------|
| `Name` | Unique identifier for each antibody | Yes | `100F4`, `K77-1A06` |
| `VH_nuc` | Heavy chain nucleotide sequence | Yes | `ATG...` |
| `VH_AA` | Heavy chain amino acid sequence | Yes | `QVQL...` |
| `VL_nuc` | Light chain nucleotide sequence | Yes | `ATG...` |
| `VL_AA` | Light chain amino acid sequence | Yes | `DIQMT...` |
| `Heavy_V_gene` | Heavy chain V gene annotation | No | `IGHV4-61*03` |
| `Heavy_J_gene` | Heavy chain J gene annotation | No | `IGHJ4*02` |
| `Heavy_D_gene` | Heavy chain D gene annotation | No | `IGHD4-17*01` |
| `Light_V_gene` | Light chain V gene annotation | No | `IGLV1-40*01` |
| `Light_J_gene` | Light chain J gene annotation | No | `IGLJ1*01` |
| `Specificity` | Antibody specificity | No | `HA:Unk`, `Group 1` |

### **Data Quality Requirements**
- **Complete sequences**: Both heavy and light chain sequences should be present
- **Valid characters**: Nucleotide sequences should contain only A, T, G, C
- **No stop codons**: Amino acid sequences should not contain `*` characters
- **Proper length**: Sequences should be of appropriate length for antibody chains
- **Unique names**: Each antibody should have a unique identifier

### **Example Input File**

NOTE: It is important to skip the first row when making the table. This is due to the special format of our table, where the first row is the title of the table

| INPUT TABLE |
| Name | VH_nuc | VH_AA | VL_nuc | VL_AA | Heavy_V_gene | Heavy_J_gene | Heavy_D_gene | Light_V_gene | Light_J_gene | Specificity |
|------|--------|-------|--------|-------|--------------|--------------|--------------|--------------|--------------|-------------|
| 100F4 | CAG... | QVQL... | CAG... | DIQMT... | IGHV4-61*03 | IGHJ4*02 | IGHD4-17*01 | IGLV1-40*01 | IGLJ1*01 | HA:Unk |
| K77-1A06 | ATG... | QVQL... | ATG... | DIQMT... | IGHV1-69 | IGHJ1 | IGHD3-9 | IGLV1 | IGLJ1 | Group 1 |

## 🔬 **Pipeline Workflow**

### **1. Data Filtering & Extraction**

**Command:**
```bash
python HA_screen/script/extract.py -i HA_screen/data/TableS1.csv -v IGHV1-69 IGHV6-1 IGHV1-18 -d IGHD3-9 -g ${pyir_db}/Ig/human -o ui_results/filtered_output.csv

# With optional clonotype filtering disabled (3x more sequences):
python HA_screen/script/extract.py -i HA_screen/data/TableS1.csv --skip-clonotype-filter -v IGHV1-69 IGHV6-1 IGHV1-18 -d IGHD3-9 -g ${pyir_db}/Ig/human -o ui_results/filtered_output.csv
```

**Parameters:**
- `-i`: Input table with antibody data
- `-v`: V gene families to filter out (e.g., IGHV1-69)
- `-d`: D gene families to filter out (e.g., IGHD3-9)
- `-g`: PyIR germline database path (auto-detected in UI)
- `-o`: Output filtered CSV file
- `--skip-clonotype-filter`: Optional flag to retain 3x more sequences

**Process:**
- Downloads and filters antibody data from research papers
- Removes unpaired, incomplete antibodies
- Uses PyIR package for annotation and germline retrieval
- Applies Kabat numbering for sequence completion
- Result: 302 sequences (with filter) or 958 sequences (without filter)

### **2. Sequence Segments Iteration**

**Command:**
```bash
python HA_screen/script/iteration.py -i ui_results/filtered_output.csv -p 2000000 -n HA_screen/result/random_neg.csv -o ui_results/iteration_output.fa
```

**Parameters:**
- `-i`: Filtered table as input (result from Step 1)
- `-p`: Total number of random sequences to generate (pool size)
- `-o`: Output FASTA file containing truncated sequences
- `-n`: Negative control sequence list

**Process:**
- Generates random codon assignments for amino acid sequences
- Creates millions of DNA sequence variants
- Truncates sequences into manageable segments (99 bp each)
- Result: 15.9M+ DNA sequences in FASTA format (1.7GB file)

### **3. CD-HIT Clustering**

**Command:**
```bash
bash HA_screen/script/cd-hit.sh
```

**Process:**
- Groups sequences based on similarity using CD-HIT algorithm
- Input: FASTA file from Step 2
- Output: Clustered sequences in `cdhit/` directory
- Note: Computationally intensive step

### **4. CD-HIT Result Selection**

**Command:**
```bash
python HA_screen/script/cdhit_result.py -i ui_results/iteration_output.fa -n HA_screen/result/random_neg.csv -gs 25 -ng 12 -nn 2
```

**Parameters:**
- `-gs`: Group size (25)
- `-ng`: Number of groups (12)
- `-nn`: Number of negative controls (2)

### **5. ChunkByOverlap Processing**

**Command:**
```bash
python HA_screen/script/ChunkByOverlap.py
```

**Process:**
- Further processes clustered sequences
- Optimizes overlap regions for synthesis

### **6. Final Pool Assembly**

**Result:**
- Production-ready DNA sequences for antibody library synthesis
- Optimized for high-throughput screening
- Compatible with standard synthesis platforms

## ✅ **Successful Setup Output Examples**

### **Python Setup Script (`python setup_cross_platform.py`)**

```bash
🚀 oPool Design Pipeline - Cross-Platform Setup
==================================================
🖥️  Platform: Darwin 23.4.0
🐍 Python: 3.11.7 (main, Dec 15 2023, 12:09:04) [Clang 14.0.6 ]

✅ Conda found: conda 24.11.3
🔧 Creating conda environment 'oPool'...

## Package Plan ##
  environment location: /Users/username/envs/oPool
  
  added / updated specs:
    - python=3.9

The following NEW packages will be INSTALLED:
  ca-certificates    conda-forge/osx-64::ca-certificates-2025.7.15
  python             conda-forge/osx-64::python-3.9.23
  pip                conda-forge/noarch::pip-25.2
  # ... additional packages

✅ Environment 'oPool' created successfully
📦 Installing dependencies...

Installing pip dependencies:
Collecting Flask==2.3.3
  Using cached flask-2.3.3-py3-none-any.whl.metadata (3.6 kB)
Collecting crowelab-pyir
  Using cached crowelab_pyir-1.5.0-py3-none-any.whl
# ... dependency installation

✅ Dependencies installed successfully
🔬 Setting up PyIR...
✅ PyIR installed successfully
🔍 Setting up germline database...
⚠️  Note: IMGT.org downloads can be slow - this may take 20+ minutes
✅ Germline database setup completed successfully
🔍 Creating directories...
🔍 Exists: oPool_design/uploads
🔍 Exists: oPool_design/results  
🔍 Exists: oPool_design/logs

🎉 Setup completed!

📋 Next steps:
1. Activate environment: conda activate oPool
2. Navigate to oPool_design: cd oPool_design
3. Launch UI: python launch_ui.py

💡 If you encounter issues, check the troubleshooting guide in CROSS_PLATFORM_SETUP.md
```

### **Python Setup Script (`./setup_cross_platform.py`)**

```bash
🚀 oPool Design Pipeline - Cross-Platform Setup
==================================================
🖥️  Platform: Darwin (macOS)
🐍 Python: 3.11.7

✅ Conda found: conda 24.11.3
🔧 Creating conda environment 'oPool'...
📦 Installing dependencies from environment.yml...
🔬 Setting up PyIR germline database...
🔍 Creating project directories...

✅ Setup completed successfully!

🎯 Quick Start:
   conda activate oPool
   cd oPool_design  
   python launch_ui.py

🔗 Access UI at: http://127.0.0.1:5001
```

### **Environment Activation & UI Launch**

```bash
# Activate the conda environment
$ conda activate oPool
(oPool) $ cd oPool_design

# Launch the web UI
(oPool) $ python launch_ui.py
🔬 oPool Design Pipeline Web UI Launcher
==================================================
📁 Current working directory: /path/to/oPool_design
✅ Found UI directory: /path/to/oPool_design/ui
✅ Flask app imported successfully
🌐 Starting web server...
🎉 Web UI launched successfully!
🔗 Access at: http://127.0.0.1:5001

# UI is now running and accessible in your browser
```

### **What Each Setup Component Does**

| Component | Purpose | Expected Output |
|-----------|---------|----------------|
| **Platform Detection** | Identifies OS (macOS/Linux) | `🖥️ Platform: Darwin 23.4.0` |
| **Conda Check** | Verifies conda installation | `✅ Conda found: conda 24.11.3` |
| **Environment Creation** | Creates isolated Python env | `✅ Environment 'oPool' created successfully` |
| **Dependency Installation** | Installs scientific packages | `📦 Installing blast, cd-hit, hmmer, pandas...` |
| **PyIR Setup** | Downloads germline database | `✅ Germline database setup completed` |
| **Directory Creation** | Creates uploads/results folders | `🔍 Created: oPool_design/uploads` |

### **Troubleshooting Common Setup Issues**

#### **If PyIR Setup Times Out:**
```bash
# The setup will automatically detect partial downloads
⏰ PyIR setup timed out (1 hour)
🔄 Checking for partial setup...
✅ Partial setup detected - proceeding with available data
```

#### **If Environment Already Exists:**
```bash
⚠️  Environment 'oPool' already exists
Do you want to recreate it? (y/N): y
🗑️  Removing existing environment 'oPool'...
✅ Environment recreated successfully
```

#### **Successful Dependency Verification:**
```bash
🔍 Found germline data for: human, mouse, rat
✅ All required dependencies installed and verified
🎯 Pipeline ready for use!
```

## 🧪 **Verified Test Results**

### **Extract**
| Test Dataset | Input Sequences | With Clonotype Filter | Without Clonotype Filter |
|--------------|-----------------|----------------------|--------------------------|
| TableS1 (Full) | 5,563 sequences | 302 sequences | **958 sequences** |
| Processing Time | - | ~1 minute | ~1 minute |

### **Iteration Step**
| Input Dataset | Output Size | Sequences Generated | Processing Time |
|---------------|-------------|-------------------|-----------------|
| 302 sequences (filtered) | 44 MB FASTA | ~800K sequences | 2 seconds |
| **958 sequences (full)** | **1.7 GB FASTA** | **15,978,330 sequences** | **43 seconds** |

### **Real Pipeline Outputs**

#### **Extract Step Output (TableS1 without clonotype filtering):**
```bash
✅ No clonotype filtering - 958 sequences retained
🧬 Attempting Kabat numbering and sequence completion...
✅ VH Kabat numbering and completion successful
✅ VL Kabat numbering and completion successful
💾 Saving results to ui_results/TableS1_no_clonotype_filter.csv
```

#### **Iteration Step Output (Full Dataset):**
```bash
Tue Sep  2 08:57:12 2025 Loading the data
Tue Sep  2 08:57:12 2025 Star the sampling...
Processing... ████████████████████████████████████████ 100% 0:00:00
Tue Sep  2 08:57:55 2025 Saving the result...
Tue Sep  2 08:58:23 2025 Done, the result is at ui_results/TableS1_full_iteration_no_clonotype.fa

# Final output: 1.7GB FASTA file with 15,978,330 DNA sequences
```

#### **Sample DNA Sequences Generated:**
```fasta
>315-13-1B02:0-0
GACAUCCAAAUGACCCAGAGUCCGUCGAGUCUUUCCGCAUCGGUAGGCGAUCGUGUAACUAUUACGUGUCGUGCGCAGCCAAACAAUUAGCCGCUAUCUU
>315-13-1B02:1-0
AAUUGGUAUCAGCAGAAAGCGGGGAAGGCGCCCACCUUACUUAUUUACGAUGCAUCGCGCCUUCAGUCCGGUGUCCCUUCUCGUUUCUCAGGGAGUGGU
>315-13-1B02:2-0
UCCGGGACAGAAUUUACAUUAACUAUUUCUAGUCUUCAGCGUGAGGAUUUCGCUACGUACUAUUGUCAGCAAUCAGACUCCAUUCCAGCUUUAACAUUC
```


### **Optional Clonotype Filtering**
```python
# Command line usage:
python HA_screen/script/extract.py -i data.csv --skip-clonotype-filter  # 3x more sequences
python HA_screen/script/extract.py -i data.csv                          # Original filtering

# Web UI: Toggle switch for "Skip Clonotype Filtering"
```

