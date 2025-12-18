import pandas as pd
import numpy as np
import os
import sys
import argparse

# Function to read raw juicer dump files
def load_juicer_dump(filename, val_col_name):
    if not os.path.exists(filename):
        print(f"Error: File not found {filename}")
        sys.exit(1)
        
    try:
        # Read data
        df = pd.read_csv(filename, sep=r'\s+', header=None, comment='#', on_bad_lines='skip')
        
        if df.shape[1] >= 5:
             df = df.iloc[:, [1, 3, 4]]
             df.columns = ['pos1', 'pos2', val_col_name]
        elif df.shape[1] == 3:
             df.columns = ['pos1', 'pos2', val_col_name]
        else:
             print(f"Error: Unknown file format {filename} (No. of columns: {df.shape[1]})")
             sys.exit(1)

        df['pos1'] = pd.to_numeric(df['pos1'], errors='coerce').fillna(0).astype(int)
        df['pos2'] = pd.to_numeric(df['pos2'], errors='coerce').fillna(0).astype(int)
        df[val_col_name] = pd.to_numeric(df[val_col_name], errors='coerce')
        
        df.dropna(inplace=True)
        return df

    except Exception as e:
        print(f"Exception found while reading {filename}: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Process Juicer dump files (Observed & OE) for DeepLoop.")
    
    parser.add_argument("--obs", required=True, help="Path to 'dump observed' file")
    parser.add_argument("--oe", required=True, help="Path to 'dump oe' file")
    parser.add_argument("--chrom", required=True, help="Chromosome name")
    parser.add_argument("--res", type=int, required=True, help="Resolution in BP")
    parser.add_argument("--out", required=True, help="Output TSV filename")
    parser.add_argument("--anchor-dir", default="anchors", help="Directory to save BED reference")

    args = parser.parse_args()

    FILE_OBS = args.obs
    FILE_OE = args.oe
    CHROM_NAME = args.chrom
    RES = args.res
    OUTPUT_FILE = args.out
    ANCHOR_DIR = args.anchor_dir

    print(f"Processing juicer dump output files for {CHROM_NAME} at {RES} resolution...")

    # Load juicer data
    df_obs = load_juicer_dump(FILE_OBS, 'observed')
    df_oe = load_juicer_dump(FILE_OE, 'oe')

    # Merge obs & oe
    merged = pd.merge(df_obs, df_oe, on=['pos1', 'pos2'], how='inner')
    print(f"Number of interactions after merge: {len(merged)}")

    # Reconstruction of expected counts
    merged['expected'] = merged['observed'] / merged['oe']
    merged.replace([np.inf, -np.inf], 0, inplace=True)
    merged.fillna(0, inplace=True)

    # Convert genomic coordinates to bin indices
    merged['anchor1'] = (merged['pos1'] / RES).astype(int)
    merged['anchor2'] = (merged['pos2'] / RES).astype(int)

    # Create final DF & sort values by anchors
    final_df = merged[['anchor1', 'anchor2', 'observed', 'expected']].copy()
    final_df.sort_values(by=['anchor1', 'anchor2'], inplace=True)

    print(f"Saving processed interaction file to: {OUTPUT_FILE}")
    final_df.to_csv(OUTPUT_FILE, sep='\t', index=False)

    # Generating BED reference
    print(f"Generating BED reference file in: {ANCHOR_DIR}...")
    if not os.path.exists(ANCHOR_DIR):
        try:
            os.makedirs(ANCHOR_DIR)
        except OSError as e:
            print(f"Error creating directory {ANCHOR_DIR}: {e}")
            sys.exit(1)

    # Determine max index
    if not final_df.empty:
        max_idx = final_df[['anchor1', 'anchor2']].max().max()
    else:
        max_idx = 0
    
    bed_file_path = os.path.join(ANCHOR_DIR, f"{CHROM_NAME}.bed")

    # BUGFIX: Increased padding from 100 to 200 to avoid error "index out of bounds"
    PADDING_BINS = 200

    try:
        with open(bed_file_path, 'w') as f:
            for i in range(int(max_idx) + PADDING_BINS):
                s = i * RES
                e = s + RES
                f.write(f"{CHROM_NAME}\t{s}\t{e}\t{i}\n")
        print(f"BED file created: {bed_file_path} (Padding: {PADDING_BINS} bins)")
    except IOError as e:
        print(f"Error writing BED file: {e}")

    print("Data ready for DeepLoop.")

if __name__ == "__main__":
    main()