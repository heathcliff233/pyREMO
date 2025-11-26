#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile

from famr_python.topology import read_famr_comm, generate_hdd, parse_pdb_sequence
from famr_python.reconstruct import run_reconstruction
from pyremo.resources import resolve_data_file

try:
    from tqdm import tqdm
except ImportError:
    # Simple fallback for tqdm
    def tqdm(iterable, desc="", total=None, unit=""):
        if total is None and hasattr(iterable, "__len__"):
            total = len(iterable)
        print(desc)
        for i, item in enumerate(iterable):
            yield item
            if i % 10 == 0 and total:
                print(f"Processed {i+1}/{total}", end='\r')
        print("")

def process_file(pdb_file, output_file, comm_data, bbdat_file, overwrite):
    if os.path.exists(output_file) and not overwrite:
        return False # Skipped

    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # 1. Parse Sequence
    try:
        seq, _ = parse_pdb_sequence(pdb_file)
        if not seq:
            # print(f"Warning: No CA atoms in {pdb_file}. Skipping.")
            return False
    except Exception as e:
        print(f"Error parsing {pdb_file}: {e}")
        return False

    # 2. Generate Hdd in temp file
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp_hdd:
        hdd_path = tmp_hdd.name
    
    try:
        generate_hdd(seq, comm_data, hdd_path)
        
        # 3. Run Reconstruction
        # Redirect stdout to devnull to keep progress bar clean
        original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        
        try:
            run_reconstruction(pdb_file, hdd_path, bbdat_file, output_file)
        except Exception as e:
            sys.stdout.close()
            sys.stdout = original_stdout
            print(f"Error reconstructing {pdb_file}: {e}")
            return False
        finally:
            if sys.stdout != original_stdout:
                sys.stdout.close()
                sys.stdout = original_stdout
            
    finally:
        if os.path.exists(hdd_path):
            os.remove(hdd_path)
            
    return True

def main():
    parser = argparse.ArgumentParser(description='REMO: Reconstruct Atomic MOdel from C-alpha trace.')
    parser.add_argument('-i', '--input', required=True, help='Input PDB file or directory.')
    parser.add_argument('-o', '--output', required=True, help='Output PDB file or directory.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output files.')
    
    args = parser.parse_args()
    
    # Resolve data assets (repo checkout or installed package)
    try:
        famr_comm = resolve_data_file('FAMRcomm')
        bbdat_file = resolve_data_file('BBdat')
    except FileNotFoundError as exc:
        print(f"Error: {exc}")
        sys.exit(1)

    print("Loading parameters...")
    comm_data = read_famr_comm(famr_comm)
    
    files_to_process = []
    
    if os.path.isdir(args.input):
        # Directory mode
        if os.path.exists(args.output) and not os.path.isdir(args.output):
            print(f"Error: Input is a directory but output {args.output} is a file.")
            sys.exit(1)
            
        input_dir = args.input
        output_dir = args.output
        
        for root, dirs, files in os.walk(input_dir):
            for file in files:
                if file.endswith('.pdb') or file.endswith('.ent'): 
                    input_path = os.path.join(root, file)
                    rel_path = os.path.relpath(input_path, input_dir)
                    output_path = os.path.join(output_dir, rel_path)
                    files_to_process.append((input_path, output_path))
                    
    elif os.path.isfile(args.input):
        # File mode
        input_path = args.input
        output_path = args.output
        
        if os.path.isdir(output_path):
            filename = os.path.basename(input_path)
            output_path = os.path.join(output_path, filename)
            
        files_to_process.append((input_path, output_path))
    else:
        print(f"Error: Input {args.input} not found.")
        sys.exit(1)
        
    print(f"Found {len(files_to_process)} files to process.")
    
    processed_count = 0
    skipped_count = 0
    
    for in_file, out_file in tqdm(files_to_process, desc="Processing", unit="file"):
        success = process_file(in_file, out_file, comm_data, bbdat_file, args.overwrite)
        if success:
            processed_count += 1
        else:
            skipped_count += 1
            
    print(f"\nDone. Processed: {processed_count}, Skipped/Error: {skipped_count}")

if __name__ == "__main__":
    main()
