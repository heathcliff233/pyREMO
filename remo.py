#!/usr/bin/env python3
import argparse
import os
import sys

from pyremo.api import load_parameters, reconstruct_file

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

def main():
    parser = argparse.ArgumentParser(description='REMO: Reconstruct Atomic MOdel from C-alpha trace.')
    parser.add_argument('-i', '--input', required=True, help='Input PDB file or directory.')
    parser.add_argument('-o', '--output', required=True, help='Output PDB file or directory.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output files.')
    
    args = parser.parse_args()
    
    print("Loading parameters...")
    # Resolve data assets (repo checkout or installed package)
    try:
        load_parameters()
    except FileNotFoundError as exc:
        print(f"Error: {exc}")
        sys.exit(1)
    
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
        success = reconstruct_file(in_file, out_file, args.overwrite)
        if success:
            processed_count += 1
        else:
            skipped_count += 1
            
    print(f"\nDone. Processed: {processed_count}, Skipped/Error: {skipped_count}")

if __name__ == "__main__":
    main()
