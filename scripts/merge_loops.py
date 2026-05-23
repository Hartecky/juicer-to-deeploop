import pandas as pd
import numpy as np
import argparse
import os
from scipy.spatial import cKDTree

def load_bedpe(filepath):
    try:
        df = pd.read_csv(filepath, sep='\t', comment='#')
        # Standard BEDPE columns
        if df.shape[1] >= 6:
            df.columns.values[:6] = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']

        # Centroids for distance calculation
        df['cx'] = (df['x1'] + df['x2']) / 2
        df['cy'] = (df['y1'] + df['y2']) / 2
        return df
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return pd.DataFrame()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--files", nargs='+', required=True, help="List of BEDPE files to merge")
    parser.add_argument("--resolutions", nargs='+', type=int, required=True, help="List of resolutions corresponding to files (e.g. 2000 5000)")
    parser.add_argument("--out", required=True, help="Output merged file")
    parser.add_argument("--tolerance", type=int, default=20000, help="Overlap tolerance in bp (Default: 20000)")
    
    args = parser.parse_args()
    
    if len(args.files) != len(args.resolutions):
        print("Error: Number of files must match number of resolutions.")
        return

    # 1. Organize data: [(2000, df_2k), (5000, df_5k), ...]
    data_list = []
    for filepath, res in zip(args.files, args.resolutions):
        print(f"Loading {res}bp: {filepath}")
        df = load_bedpe(filepath)
        if not df.empty:
            df['source_res'] = res
            data_list.append((res, df))
    
    # Sort ascending by resolution (2000 has priority over 10000)
    data_list.sort(key=lambda x: x[0])

    if not data_list:
        print("No valid loops loaded.")
        return

    # 2. Merging Algorithm
    # Base = highest resolution
    final_loops = data_list[0][1].copy()
    print(f"Base loops (from {data_list[0][0]}bp): {len(final_loops)}")

    # Iterate over remaining (lower) resolutions
    for res, candidate_df in data_list[1:]:
        if final_loops.empty:
            final_loops = pd.concat([final_loops, candidate_df])
            continue

        # Build tree from current base
        tree = cKDTree(final_loops[['cx', 'cy']].values)
        
        # Check candidates
        dists, _ = tree.query(candidate_df[['cx', 'cy']].values, k=1, distance_upper_bound=args.tolerance)
        
        # If dist == inf, it means no neighbor within range -> Unique loop
        is_unique = (dists == float('inf'))
        unique_loops = candidate_df[is_unique]
        
        if not unique_loops.empty:
            print(f"Adding {len(unique_loops)} unique loops from {res}bp")
            final_loops = pd.concat([final_loops, unique_loops], ignore_index=True)
        else:
            print(f"No unique loops in {res}bp")

    # 3. Write output
    cols = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'name', 'score', 'color', 'source_res']
    # Match available columns
    out_cols = [c for c in cols if c in final_loops.columns]
    
    final_loops = final_loops.sort_values(by=['chr1', 'x1', 'x2'])
    final_loops.to_csv(args.out, sep='\t', index=False)
    print(f"Total merged loops: {len(final_loops)}")

if __name__ == "__main__":
    main()