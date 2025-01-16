#!/usr/bin/env python3
import subprocess
import os
from datetime import datetime
import re

def convert_to_tb(size_str):
    """Convert size string (like '1.2T' or '750G') to TB float"""
    size = float(size_str[:-1])
    unit = size_str[-1].upper()
    if unit == 'T':
        return size
    elif unit == 'G':
        return size / 1024
    elif unit == 'P':
        return size * 1024
    else:
        return 0  # for smaller units, return 0 TB

def analyze_directory(path):
    """
    Analyzes a directory to get file count and total size.
    Returns a tuple of (file_count, directory_size)
    """
    try:
        # Change to the directory
        os.chdir(path)
        
        # Count files
        find_cmd = ["find", "./", "-type", "f"]
        wc_cmd = ["wc", "-l"]
        
        # Pipeline the commands
        find_process = subprocess.Popen(find_cmd, stdout=subprocess.PIPE)
        wc_process = subprocess.Popen(wc_cmd, stdin=find_process.stdout, stdout=subprocess.PIPE)
        find_process.stdout.close()
        
        file_count = int(wc_process.communicate()[0].decode().strip())
        
        # Get directory size
        du_cmd = ["du", "-sh"]
        du_output = subprocess.check_output(du_cmd).decode().strip()
        directory_size = du_output.split()[0]
        
        return file_count, directory_size
    
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        return None, None
    except Exception as e:
        print(f"Error processing directory {path}: {e}")
        return None, None

def generate_latex_table(results):
    """Generate LaTeX table from results"""
    latex_output = []
    latex_output.append(r"\begin{deluxetable*}{lccl}[b]")
    latex_output.append(r"\tablecaption{DR1 Directory Structure and Data Volume}")
    latex_output.append(r"\label{tab:directory_structure}")
    latex_output.append(r"\tablehead{")
    latex_output.append(r"\colhead{} & \colhead{Size} & \colhead{No. of} & \colhead{} \\")
    latex_output.append(r"\colhead{Directory} & \colhead{(TB)} & \colhead{Files} & \colhead{Description}")
    latex_output.append(r"}")
    latex_output.append(r"\startdata")
    
    # Add data rows
    for path, (count, size) in results.items():
        # Calculate indent level based on number of path components
        indent_level = path.count('/') - 3  # Adjust this based on your base path
        indent = r"\hspace{" + f"{0.4 * indent_level:.1f}" + r"cm}"
        
        # Format the directory name
        dir_name = path.split('/')[-1] + '/'
        dir_latex = indent + r"\texttt{" + dir_name + r"}"
        
        # Format size in TB
        size_tb = convert_to_tb(size)
        size_str = f"{size_tb:.3f}" if size_tb > 0 else "XX"
        
        # Format file count with commas
        count_str = f"{count:,}" if count is not None else "XX"
        
        # Add description placeholder
        description = "Description needed"
        
        latex_output.append(f"{dir_latex:<40} & {size_str} & {count_str} & {description} \\\\")
    
    latex_output.append(r"\enddata")
    latex_output.append(r"\tablecomments{Please see the DESI data model documentation at \url{https://desidatamodel.readthedocs.io} for more details.}")
    latex_output.append(r"\end{deluxetable*}")
    
    return "\n".join(latex_output)

def main():
    # List of directories to analyze
    base_dir = "/global/cfs/cdirs/desi/public/dr1"
    directories = [
        "spectro/data",
        "spectro/desi_spectro_calib",
        "spectro/redux/iron",
        "spectro/redux/iron/exposures",
        "spectro/redux/iron/healpix",
        "spectro/redux/iron/tiles/cumulative",
        "spectro/redux/iron/tiles/pernight",
        "spectro/redux/iron/zcatalog/v1",
        "spectro/redux/guadalupe",
        "spectro/templates",
        "survey",
        "target/catalogs",
        "target/fiberassign",
        "vac/dr1"
    ]
    
    # Store results
    results = {}
    
    # Print stdout report
    print(f"\nDirectory Analysis Report - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 80)
    print(f"{'Directory':<60} {'Files':<10} {'Size':<10}")
    print("-" * 80)
    
    for dir_path in directories:
        full_path = os.path.join(base_dir, dir_path)
        file_count, dir_size = analyze_directory(full_path)
        results[dir_path] = (file_count, dir_size)
        
        if file_count is not None and dir_size is not None:
            print(f"{full_path:<60} {file_count:<10} {dir_size:<10}")
        else:
            print(f"{full_path:<60} {'ERROR':<10} {'ERROR':<10}")
    
    print("-" * 80)
    
    # Generate and print LaTeX table to screen
    latex_table = generate_latex_table(results)
    print("\nLaTeX table content:")
    print("-" * 80)
    print(latex_table)
    print("-" * 80)

if __name__ == "__main__":
    main()