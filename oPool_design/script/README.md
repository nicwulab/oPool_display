# UI Scripts Directory

This directory contains UI-specific modified versions of the original pipeline scripts.

## Files:

- **`cdhit_result_modified.py`**: Modified version of cdhit_result.py that outputs to `ui_results/` instead of `result/`
- **`cdhit_result_backup.py`**: Backup of the original cdhit_result.py before modifications
- **`Overlap_check_modified.py`**: Modified version of Overlap_check.py that works with the new directory structure
- **`iteration_backup.py`**: Backup of iteration.py before multiprocessing fixes

## Usage:

These scripts are used for UI testing and development. They maintain compatibility with the new organized directory structure while preserving the original scripts in `HA_screen/script/`.

## Original Scripts:

The original, unmodified scripts are preserved in `HA_screen/script/` for reference and traditional command-line usage.
