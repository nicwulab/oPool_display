from flask import Flask, render_template, request, jsonify, send_file, flash, redirect, url_for, Response
import os
import subprocess
import pandas as pd
import tempfile
import shutil
from werkzeug.utils import secure_filename
import json
import threading
import time
import signal
import psutil
from pathlib import Path
import uuid

# Import configuration first
from config import get_config, get_germline_path
config = get_config()

# Create Flask app
app = Flask(__name__)
app.secret_key = config.SECRET_KEY
app.config['MAX_CONTENT_LENGTH'] = config.MAX_CONTENT_LENGTH

# Configuration
UPLOAD_FOLDER = str(config.UPLOAD_FOLDER)

# Global variables to track running processes
running_processes = {}
process_outputs = {}  # Store process outputs for streaming
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
        
        # Start the process and track it
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=parent_dir)
        running_processes[step_name] = process
        
        # Wait for the process to complete
        stdout, stderr = process.communicate()
        
        # Remove from running processes
        if step_name in running_processes:
            del running_processes[step_name]
        
        # Scan for output files
        files_created = scan_output_files(step_name, parent_dir)
        
        if process.returncode == 0:
            return {
                "success": True,
                "output": stdout,
                "step": step_name,
                "files_created": files_created
            }
        else:
            return {
                "success": False,
                "error": stderr,
                "step": step_name,
                "files_created": files_created
            }
    except Exception as e:
        # Remove from running processes if there was an error
        if step_name in running_processes:
            del running_processes[step_name]
        return {
            "success": False,
            "error": str(e),
            "step": step_name,
            "files_created": []
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
    """Run the extract pipeline step"""
    try:
        data = request.json
        input_file = data['input_file']
        v_list = data.get('v_list', ['IGHV1-69', 'IGHV6-1', 'IGHV1-18'])
        d_list = data.get('d_list', ['IGHD3-9'])
        germline_path = data.get('germline_path', '')
        skip_clonotype_filter = data.get('skip_clonotype_filter', False)
        
        # Create step1 output directory
        os.makedirs('ui_results/step1', exist_ok=True)
        
        # Construct output filename
        base_name = os.path.splitext(input_file)[0]
        output_file = f"ui_results/step1/{base_name}_output.csv"
        
        # Run extract command
        cmd_parts = [
            "python script/extract.py",
            f"-i uploads/{input_file}",
            f"-o {output_file}"
        ]
        
        if v_list:
            cmd_parts.append(f"-v {' '.join(v_list)}")
        if d_list:
            cmd_parts.append(f"-d {' '.join(d_list)}")
        if germline_path:
            cmd_parts.append(f"-g {germline_path}")
        if skip_clonotype_filter:
            cmd_parts.append("--skip-clonotype")
            
        cmd = " ".join(cmd_parts)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=os.getcwd())
        
        if result.returncode == 0:
            # Try to read and preview the output
            if os.path.exists(output_file):
                df = pd.read_csv(output_file)
                return jsonify({
                    'success': True,
                    'step': 'Extract',
                    'message': 'Extract completed successfully',
                    'output_file': output_file,
                    'total_rows': len(df),
                    'preview': df.head(10).fillna('').to_dict('records'),
                    'files_created': scan_output_files('Extract', os.getcwd())
                })
            else:
                return jsonify({
                    'success': True,
                    'step': 'Extract',
                    'message': 'Extract completed but output file not found',
                    'files_created': scan_output_files('Extract', os.getcwd())
                })
        else:
            return jsonify({
                'success': False,
                'step': 'Extract',
                'error': result.stderr or 'Unknown error occurred',
                'files_created': scan_output_files('Extract', os.getcwd())
            })
            
    except Exception as e:
        return jsonify({
            'success': False,
            'step': 'Extract',
            'error': str(e),
            'files_created': []
        })

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
    """Run the CD-HIT clustering step and return process ID for streaming"""
    try:
        data = request.get_json()
        
        # Get parameters
        min_threshold = data.get('min_threshold', 0.6)
        max_threshold = data.get('max_threshold', 0.85)
        increment = data.get('increment', 0.05)
        
        # Generate unique process ID
        process_id = str(uuid.uuid4())
        
        # Run cd-hit script with parameters
        cmd = f"bash script/cd-hit.sh {min_threshold} {max_threshold} {increment}"
        
        # Start the process in the background
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        process = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=parent_dir,
            bufsize=0,  # Unbuffered
            universal_newlines=True
        )
        
        # Store process information
        running_processes[process_id] = process
        process_outputs[process_id] = []
        
        # Start a thread to read output
        def read_output():
            try:
                for line in iter(process.stdout.readline, ''):
                    if line:
                        # Keep line breaks for proper display
                        process_outputs[process_id].append(line.rstrip('\n\r'))
                process.stdout.close()
                process.wait()
                # Mark process as completed
                if process_id in running_processes:
                    del running_processes[process_id]
            except Exception as e:
                # Handle any errors in reading output
                process_outputs[process_id].append(f"Error reading output: {str(e)}")
                if process_id in running_processes:
                    del running_processes[process_id]
        
        thread = threading.Thread(target=read_output)
        thread.daemon = True
        thread.start()
        
        return jsonify({
            'success': True,
            'process_id': process_id,
            'message': 'CD-HIT process started successfully'
        })
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })

@app.route('/stream_output/<process_id>')
def stream_output(process_id):
    """Stream real-time output from a running process"""
    def generate():
        last_index = 0
        while True:
            if process_id in process_outputs:
                current_outputs = process_outputs[process_id]
                if len(current_outputs) > last_index:
                    # Send new lines with the format the frontend expects
                    for i in range(last_index, len(current_outputs)):
                        content = current_outputs[i] + '\n'
                        yield f"data: {json.dumps({'type': 'output', 'content': content})}\n\n"
                    last_index = len(current_outputs)
                
                # Check if process is still running
                if process_id not in running_processes:
                    # Process completed, send final status
                    yield f"data: {json.dumps({'type': 'done', 'return_code': 0})}\n\n"
                    break
            else:
                # Process not found
                yield f"data: {json.dumps({'type': 'error', 'content': 'Process not found'})}\n\n"
                break
            
            time.sleep(0.1)  # Check more frequently (every 100ms)
    
    return Response(generate(), mimetype='text/event-stream')

@app.route('/abort_process/<process_id>', methods=['POST'])
def abort_process(process_id):
    """Abort a running process"""
    try:
        if process_id in running_processes:
            process = running_processes[process_id]
            process.terminate()
            del running_processes[process_id]
            if process_id in process_outputs:
                del process_outputs[process_id]
            return jsonify({'success': True, 'message': 'Process aborted'})
        else:
            return jsonify({'success': False, 'error': 'Process not found'})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@app.route('/run_cdhit_result', methods=['POST'])
def run_cdhit_result():
    """Run the CD-HIT result selection step"""
    try:
        data = request.json
        negative_file = data.get('negative_file')
        group_size = data.get('group_size', 25)
        num_groups = data.get('num_groups', 12)
        num_negative = data.get('num_negative', 2)
        
        # Create step4 output directory
        os.makedirs('ui_results/step4', exist_ok=True)
        
        # Save group_size for Step 5 to use
        with open('ui_results/step4/group_size.txt', 'w') as f:
            f.write(str(group_size))
        
        # Input: Always use the most recent FASTA file from step2 (segmented sequences)
        step2_dir = "ui_results/step2"
        input_path = None
        if os.path.exists(step2_dir):
            fasta_files = [f for f in os.listdir(step2_dir) if f.endswith('.fa')]
            if fasta_files:
                fasta_files.sort(key=lambda x: os.path.getmtime(os.path.join(step2_dir, x)), reverse=True)
                input_path = f"ui_results/step2/{fasta_files[0]}"
        
        if not input_path:
            return jsonify({
                'success': False,
                'step': 'CD-HIT Result Selection',
                'error': 'No FASTA file found in step2. Please complete Step 2 (Iteration) first.',
                'files_created': []
            })
            
        # Negative: Use user-selected CSV file
        if not negative_file:
            return jsonify({
                'success': False,
                'step': 'CD-HIT Result Selection',
                'error': 'No negative control file specified. Please select a negative control file.',
                'files_created': []
            })
        
        # Find the negative file in appropriate directories
        negative_path = None
        for search_dir in ['ui_results/step1', 'uploads']:
            test_path = f"{search_dir}/{negative_file}"
            if os.path.exists(test_path):
                negative_path = test_path
                break
        
        if not negative_path:
            return jsonify({
                'success': False,
                'step': 'CD-HIT Result Selection',
                'error': f'Negative control file "{negative_file}" not found. Please check the file exists.',
                'files_created': []
            })
        
        cmd = f"python script/cdhit_result_modified.py -i {input_path} -n {negative_path} -gs {group_size} -ng {num_groups} -nn {num_negative}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=os.getcwd())
        
        if result.returncode == 0:
            return jsonify({
                'success': True,
                'step': 'CD-HIT Result Selection',
                'message': f'CD-HIT result selection completed successfully using {os.path.basename(input_path)} and {os.path.basename(negative_path)}',
                'files_created': scan_output_files('cdhit_result', os.getcwd()),
                'input_used': input_path,
                'negative_used': negative_path
            })
        else:
            return jsonify({
                'success': False,
                'step': 'CD-HIT Result Selection',
                'error': result.stderr or result.stdout or 'Unknown error occurred',
                'files_created': scan_output_files('cdhit_result', os.getcwd()),
                'input_attempted': input_path,
                'negative_attempted': negative_path
            })
            
    except Exception as e:
        return jsonify({
            'success': False,
            'step': 'CD-HIT Result Selection',
            'error': str(e),
            'files_created': []
        })

@app.route('/run_overlap_check', methods=['POST'])
def run_overlap_check():
    """Run the overlap check step"""
    try:
        data = request.json
        input_file = data.get('input_file', '')  # This parameter is not used by the script
        negative_file = data['negative_file']
        
        # Create step5 output directories
        os.makedirs('ui_results/step5/primers', exist_ok=True)
        os.makedirs('ui_results/step5/blast', exist_ok=True)
        os.makedirs('ui_results/step5/segs_id', exist_ok=True)
        
        # Find negative file
        negative_path = f"ui_results/step1/{negative_file}" if not negative_file.startswith('ui_results') else negative_file
        if not os.path.exists(negative_path):
            negative_path = f"uploads/{negative_file}"
        
        if not os.path.exists(negative_path):
            return jsonify({
                'success': False,
                'step': 'Overlap Check',
                'error': f'Negative control file "{negative_file}" not found',
                'files_created': []
            })
        
        # Try to read the group_size from Step 4
        group_size = 25  # default
        try:
            with open('ui_results/step4/group_size.txt', 'r') as f:
                group_size = int(f.read().strip())
        except (FileNotFoundError, ValueError):
            pass
        
        cmd = f"python script/Overlap_check_modified.py -n {negative_path} -g {group_size}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=os.getcwd())
        
        if result.returncode == 0:
            return jsonify({
                'success': True,
                'step': 'Overlap Check',
                'message': f'Overlap check completed successfully using {negative_file}',
                'files_created': scan_output_files('overlap_check', os.getcwd()),
                'negative_used': negative_path
            })
        else:
            return jsonify({
                'success': False,
                'step': 'Overlap Check',
                'error': result.stderr or 'Unknown error occurred',
                'files_created': scan_output_files('overlap_check', os.getcwd()),
                'negative_attempted': negative_path
            })
            
    except Exception as e:
        return jsonify({
            'success': False,
            'step': 'Overlap Check',
            'error': str(e),
            'files_created': []
        })

@app.route('/run_chunk_by_overlap', methods=['POST'])
def run_chunk_by_overlap():
    data = request.get_json()
    
    cmd = f"python script/ChunkByOverlap.py"
    
    result = run_pipeline_step('Chunk by Overlap', cmd)
    
    return jsonify(result)

@app.route("/get_overlap_input_files")
def get_overlap_input_files():
    """Get list of input files available for Step 5 (Overlap Check) from Step 4 output"""
    try:
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        ui_results_dir = os.path.join(parent_dir, "ui_results")
        
        input_files = []
        if os.path.exists(ui_results_dir):
            # Look for Re_assembled_*.fa files in step4 directory (output from cdhit_result_modified.py)
            step4_dir = os.path.join(ui_results_dir, "step4")
            if os.path.exists(step4_dir):
                for file in os.listdir(step4_dir):
                    if file.startswith("Re_assembled_") and file.endswith(".fa"):
                        input_files.append(file)
        
        return jsonify({
            "success": True,
            "input_files": input_files,
            "count": len(input_files)
        })
    except Exception as e:
        return jsonify({
            "success": False,
            "error": str(e),
            "input_files": [],
            "count": 0
        })

@app.route('/get_file_list')
def get_file_list():
    """Get list of available files in upload, result, and ui_results folders"""
    upload_files = []
    result_files = []
    ui_results_files = []
    
    # Scan upload folder
    for filename in os.listdir(app.config['UPLOAD_FOLDER']):
        if allowed_file(filename):
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            upload_files.append({
                'name': filename,
                'path': filepath,
                'size': os.path.getsize(filepath),
                'modified': time.ctime(os.path.getmtime(filepath)),
                'step': 'uploads'
            })
    
    # Scan result folder
    for filename in os.listdir(app.config['RESULT_FOLDER']):
        if allowed_file(filename):
            filepath = os.path.join(app.config['RESULT_FOLDER'], filename)
            result_files.append({
                'name': filename,
                'path': filepath,
                'size': os.path.getsize(filepath),
                'modified': time.ctime(os.path.getmtime(filepath)),
                'step': 'results'
            })
    
    # Scan ui_results directories
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ui_results_dir = os.path.join(parent_dir, "ui_results")
    
    if os.path.exists(ui_results_dir):
        for step_dir in os.listdir(ui_results_dir):
            step_path = os.path.join(ui_results_dir, step_dir)
            if os.path.isdir(step_path):
                for filename in os.listdir(step_path):
                    if allowed_file(filename):
                        filepath = os.path.join(step_path, filename)
                        ui_results_files.append({
                            'name': filename,
                            'path': filepath,
                            'size': os.path.getsize(filepath),
                            'modified': time.ctime(os.path.getmtime(filepath)),
                            'step': step_dir
                        })
    
    return jsonify({
        'upload_files': upload_files,
        'result_files': result_files,
        'ui_results_files': ui_results_files
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
