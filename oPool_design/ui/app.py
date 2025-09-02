from flask import Flask, render_template, request, jsonify, send_file, flash, redirect, url_for
import os
import subprocess
import pandas as pd
import tempfile
import shutil
from werkzeug.utils import secure_filename
import json
import threading
import time
from pathlib import Path

# Import configuration first
from config import get_config, get_germline_path
config = get_config()

# Create Flask app
app = Flask(__name__)
app.secret_key = config.SECRET_KEY
app.config['MAX_CONTENT_LENGTH'] = config.MAX_CONTENT_LENGTH

# Configuration
UPLOAD_FOLDER = str(config.UPLOAD_FOLDER)
RESULT_FOLDER = str(config.RESULT_FOLDER)
ALLOWED_EXTENSIONS = config.ALLOWED_EXTENSIONS

# Create directories if they don't exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULT_FOLDER, exist_ok=True)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['RESULT_FOLDER'] = RESULT_FOLDER
app.config['SECRET_KEY'] = config.SECRET_KEY

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def get_file_size_readable(filepath):
    """Get human-readable file size"""
    try:
        size = os.path.getsize(filepath)
        for unit in ["B", "KB", "MB", "GB"]:
            if size < 1024.0:
                return f"{size:.1f} {unit}"
            size /= 1024.0
        return f"{size:.1f} TB"
    except OSError:
        return "Unknown"

def scan_output_files(step_name, parent_dir):
    """Scan for output files created by each pipeline step"""
    files_created = []
    
    if step_name == "Extract":
        expected_files = [
            ("ui_results/VH.fa", "Heavy chain sequences (FASTA)", "Intermediate"),
            ("ui_results/VL.fa", "Light chain sequences (FASTA)", "Intermediate"),
            ("ui_results/VH.tsv.gz", "Heavy chain PyIR annotations", "Intermediate"),
            ("ui_results/VL.tsv.gz", "Light chain PyIR annotations", "Intermediate")
        ]
        try:
            for file in os.listdir(os.path.join(parent_dir, "ui_results")):
                if file.endswith("_output.csv"):
                    expected_files.append((f"ui_results/{file}", "Filtered antibody data (Main Output)", "Primary"))
        except:
            pass
    elif step_name == "Iteration":
        expected_files = []
        try:
            for file in os.listdir(os.path.join(parent_dir, "ui_results")):
                if file.startswith("iteration_") and file.endswith(".fa"):
                    expected_files.append((f"ui_results/{file}", "Segmented sequences with random codons (Main Output)", "Primary"))
        except:
            pass
    elif step_name == "CD-HIT":
        expected_files = [("cdhit/cd-hit.log", "CD-HIT clustering log", "Log")]
        try:
            for file in os.listdir(os.path.join(parent_dir, "ui_results")):
                if file.endswith(".clstr"):
                    expected_files.append((f"ui_results/{file}", "Sequence cluster file", "Intermediate"))
        except:
            pass
    elif step_name == "Overlap Check":
        expected_files = [("ui_results/segs_id/", "Segment overlap analysis files", "Directory")]
        try:
            segs_dir = os.path.join(parent_dir, "ui_results/segs_id")
            if os.path.exists(segs_dir):
                for file in os.listdir(segs_dir):
                    if file.endswith(".csv"):
                        expected_files.append((f"ui_results/segs_id/{file}", "Overlap analysis results", "Intermediate"))
        except:
            pass
    elif step_name == "Chunk by Overlap":
        expected_files = []
        try:
            for file in os.listdir(os.path.join(parent_dir, "ui_results")):
                if "_oPool_design.csv" in file:
                    expected_files.append((f"ui_results/{file}", "Final library design (Main Output)", "Primary"))
                elif "_final_sequences.tsv" in file:
                    expected_files.append((f"ui_results/{file}", "Final sequences with metadata", "Primary"))
                elif file.endswith("_sequences.fasta"):
                    expected_files.append((f"ui_results/{file}", "Pool-specific FASTA sequences", "Primary"))
        except:
            pass
    else:
        expected_files = []
    
    for filepath, description, file_type in expected_files:
        full_path = os.path.join(parent_dir, filepath)
        if os.path.exists(full_path):
            if os.path.isdir(full_path):
                try:
                    file_count = len([f for f in os.listdir(full_path) if os.path.isfile(os.path.join(full_path, f))])
                    files_created.append({"path": filepath, "description": description, "type": file_type, "size": f"{file_count} files", "exists": True})
                except:
                    files_created.append({"path": filepath, "description": description, "type": file_type, "size": "Unknown", "exists": True})
            else:
                files_created.append({"path": filepath, "description": description, "type": file_type, "size": get_file_size_readable(full_path), "exists": True})
        else:
            files_created.append({"path": filepath, "description": description, "type": file_type, "size": "Not created", "exists": False})
    
    return files_created

def run_pipeline_step(step_name, command, output_file=None):
    """Run a pipeline step and return results"""
    try:
        # Run command from the parent directory (oPool_design) not ui/
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        result = subprocess.run(command, shell=True, capture_output=True, text=True, cwd=parent_dir)
        
        # Scan for output files
        files_created = scan_output_files(step_name, parent_dir)
        
        if result.returncode == 0:
            return {
                'success': True,
                'output': result.stdout,
                'step': step_name,
                'files_created': files_created
            }
        else:
            return {
                'success': False,
                'error': result.stderr,
                'step': step_name,
                'files_created': files_created
            }
    except Exception as e:
        return {
            'success': False,
            'error': str(e),
            'step': step_name,
            'files_created': []
        }

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400
    
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400
    
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)
        
        return jsonify({
            'success': True,
            'filename': filename,
            'filepath': filepath
        })
    
    return jsonify({'error': 'Invalid file type'}), 400

@app.route('/run_extract', methods=['POST'])
def run_extract():
    data = request.get_json()
    
    # Build command for extract.py
    input_file = data['input_file']
    # If it's just a filename, assume it's in uploads directory
    if not os.path.isabs(input_file):
        input_file = os.path.join(app.config['UPLOAD_FOLDER'], input_file)
    
    # Make sure the file exists
    if not os.path.exists(input_file):
        return jsonify({
            'success': False,
            'error': f'Input file not found: {input_file}',
            'step': 'Extract'
        })
    
    cmd = f"python script/extract.py -i {input_file}"
    
    v_list = data.get('v_list', config.DEFAULT_V_GENE_FAMILIES)
    d_list = data.get('d_list', config.DEFAULT_D_GENE_FAMILIES)
    
    if v_list:
        cmd += f" -v {' '.join(v_list)}"
    
    if d_list:
        cmd += f" -d {' '.join(d_list)}"
    
    # Get germline path - use provided path or auto-detect
    germline_path = data.get('germline_path')
    if not germline_path:
        germline_path = get_germline_path()
        if not germline_path:
            return jsonify({
                'success': False,
                'error': 'No germline database path found. Please set GERMLINE_DB_PATH environment variable or install PyIR.',
                'step': 'Extract'
            })
    
    if germline_path:
        cmd += f" -g {germline_path}"
    
    # Add clonotype filter parameter
    skip_clonotype_filter = data.get('skip_clonotype_filter', config.DEFAULT_SKIP_CLONOTYPE_FILTER)
    if skip_clonotype_filter:
        cmd += " --skip-clonotype-filter"
    
    output_file = os.path.join(app.config['RESULT_FOLDER'], 'extract_output.csv')
    cmd += f" -o {output_file}"
    
    result = run_pipeline_step('Extract', cmd, output_file)
    
    if result['success']:
        # Try to read the output file to show preview
        try:
            df = pd.read_csv(output_file)
            # Replace NaN values with None for JSON compatibility
            df_preview = df.head(10).fillna('')
            result['preview'] = df_preview.to_dict('records')
            result['total_rows'] = len(df)
        except:
            result['preview'] = []
            result['total_rows'] = 0
    
    return jsonify(result)

@app.route('/run_iteration', methods=['POST'])
def run_iteration():
    data = request.get_json()
    
    pool_size = data.get('pool_size', config.DEFAULT_POOL_SIZE)
    cmd = f"python script/iteration.py -i {data['input_file']} -p {pool_size} -n {data['negative_file']} -o {data['output_file']}"
    
    result = run_pipeline_step('Iteration', cmd, data['output_file'])
    
    if result['success']:
        # Check if output file exists and get file size
        if os.path.exists(data['output_file']):
            result['file_size'] = os.path.getsize(data['output_file'])
            result['file_size_mb'] = round(result['file_size'] / (1024*1024), 2)
    
    return jsonify(result)

@app.route('/run_cdhit', methods=['POST'])
def run_cdhit():
    data = request.get_json()
    
    # Run cd-hit script
    cmd = f"bash script/cd-hit.sh"
    
    result = run_pipeline_step('CD-HIT', cmd)
    
    return jsonify(result)

@app.route('/run_cdhit_result', methods=['POST'])
def run_cdhit_result():
    data = request.get_json()
    
    group_size = data.get('group_size', config.DEFAULT_GROUP_SIZE)
    num_groups = data.get('num_groups', config.DEFAULT_NUM_GROUPS)
    num_negative = data.get('num_negative', config.DEFAULT_NUM_NEGATIVE)
    cmd = f"python script/cdhit_result_modified.py -i {data['input_file']} -n {data['negative_file']} -gs {group_size} -ng {num_groups} -nn {num_negative}"
    
    result = run_pipeline_step('CD-HIT Result Selection', cmd)
    
    return jsonify(result)

@app.route('/run_overlap_check', methods=['POST'])
def run_overlap_check():
    data = request.get_json()
    
    cmd = f"python script/Overlap_check_modified.py -i {data['input_file']} -n {data['negative_file']}"
    
    result = run_pipeline_step('Overlap Check', cmd)
    
    return jsonify(result)

@app.route('/run_chunk_by_overlap', methods=['POST'])
def run_chunk_by_overlap():
    data = request.get_json()
    
    cmd = f"python script/ChunkByOverlap.py"
    
    result = run_pipeline_step('Chunk by Overlap', cmd)
    
    return jsonify(result)

@app.route('/get_file_list')
def get_file_list():
    """Get list of available files in upload and result folders"""
    upload_files = []
    result_files = []
    
    for filename in os.listdir(app.config['UPLOAD_FOLDER']):
        if allowed_file(filename):
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            upload_files.append({
                'name': filename,
                'path': filepath,
                'size': os.path.getsize(filepath),
                'modified': time.ctime(os.path.getmtime(filepath))
            })
    
    for filename in os.listdir(app.config['RESULT_FOLDER']):
        if allowed_file(filename):
            filepath = os.path.join(app.config['RESULT_FOLDER'], filename)
            result_files.append({
                'name': filename,
                'path': filepath,
                'size': os.path.getsize(filepath),
                'modified': time.ctime(os.path.getmtime(filepath))
            })
    
    return jsonify({
        'upload_files': upload_files,
        'result_files': result_files
    })

@app.route('/download/<filename>')
def download_file(filename):
    """Download a file from the results folder"""
    try:
        return send_file(
            os.path.join(app.config['RESULT_FOLDER'], filename),
            as_attachment=True,
            download_name=filename
        )
    except FileNotFoundError:
        return jsonify({'error': 'File not found'}), 404

@app.route('/view_file/<path:filepath>')
def view_file(filepath):
    """View any file in the pipeline directories"""
    try:
        # Build the full path - filepath already includes directory structure
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        full_path = os.path.join(parent_dir, filepath)
        
        if not os.path.exists(full_path):
            return jsonify({'error': 'File not found'}), 404
        
        filename = os.path.basename(filepath)
        
        if filename.endswith(('.csv', '.tsv')):
            # For CSV/TSV files, read with pandas
            if filename.endswith('.tsv'):
                df = pd.read_csv(full_path, sep='\t')
            else:
                df = pd.read_csv(full_path)
            
            return jsonify({
                'type': 'table',
                'data': df.head(100).fillna('').to_dict('records'),
                'columns': df.columns.tolist(),
                'total_rows': len(df),
                'file_path': filepath
            })
        
        elif filename.endswith(('.fa', '.fasta')):
            # For FASTA files, read first 100 lines
            with open(full_path, 'r') as f:
                lines = f.readlines()[:100]
            
            return jsonify({
                'type': 'fasta',
                'data': lines,
                'total_lines': len(lines),
                'file_path': filepath
            })
        
        else:
            # For other files, read as text
            with open(full_path, 'r') as f:
                content = f.read(5000)  # First 5000 characters
            
            return jsonify({
                'type': 'text',
                'data': content,
                'truncated': len(content) == 5000,
                'file_path': filepath
            })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/preview/<filename>')
def preview_file(filename):
    """Preview a file (first 100 lines)"""
    try:
        filepath = os.path.join(app.config['RESULT_FOLDER'], filename)
        
        if filename.endswith(('.csv', '.tsv')):
            # For CSV/TSV files, read with pandas
            if filename.endswith('.tsv'):
                df = pd.read_csv(filepath, sep='\t')
            else:
                df = pd.read_csv(filepath)
            
            return jsonify({
                'type': 'table',
                'data': df.head(100).fillna('').to_dict('records'),
                'columns': df.columns.tolist(),
                'total_rows': len(df)
            })
        
        elif filename.endswith(('.fa', '.fasta')):
            # For FASTA files, read first 100 lines
            with open(filepath, 'r') as f:
                lines = f.readlines()[:100]
            
            return jsonify({
                'type': 'fasta',
                'data': lines,
                'total_lines': len(lines)
            })
        
        else:
            # For other files, read as text
            with open(filepath, 'r') as f:
                content = f.read(5000)  # First 5000 characters
            
            return jsonify({
                'type': 'text',
                'data': content,
                'truncated': len(content) == 5000
            })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/get_germline_path')
def get_germline_path_route():
    """Get the detected germline database path"""
    path = get_germline_path()
    if path:
        return jsonify({
            'success': True,
            'path': path,
            'message': 'Germline database path detected automatically'
        })
    else:
        return jsonify({
            'success': False,
            'path': '',
            'message': 'No germline database path found. Please set GERMLINE_DB_PATH environment variable or install PyIR.'
        })

@app.route('/environment_info')
def environment_info():
    """Get environment setup information"""
    from config import get_environment_setup_instructions
    instructions = get_environment_setup_instructions()
    return jsonify({'instructions': instructions})

if __name__ == '__main__':
    app.run(debug=config.DEBUG, host=config.HOST, port=config.PORT)
