import pandas as pd
import numpy as np
import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="Convert DeepLoop output to BEDPE format (Percentile only).")
    
    parser.add_argument("--input", required=True, help="Path to DeepLoop output file")
    parser.add_argument("--out", required=True, help="Output BEDPE filename")
    parser.add_argument("--chrom", required=True, help="Chromosome name")
    parser.add_argument("--res", type=int, required=True, help="Resolution in BP")
    
    # Default parameters
    parser.add_argument("--min-dist", type=int, default=5, help="Min distance in bins (Default: 5)")
    parser.add_argument("--percentile", type=float, default=98.0, help="Percentile (Default: 98.0)")

    args = parser.parse_args()

    # Capture variables from argparse
    INPUT_FILE = args.input
    OUTPUT_FILE = args.out
    CHROM = args.chrom
    RES = args.res
    MIN_BIN_DIST = args.min_dist
    PERCENTILE = args.percentile

    print(f"DeepLoop -> BEDPE: {CHROM} @ {RES}bp")
    print(f"Strategy: Percentile ({PERCENTILE}%) | MinBinDist={MIN_BIN_DIST}")

    try:
        if not os.path.exists(INPUT_FILE):
            raise FileNotFoundError(f"Missing file: {INPUT_FILE}")
        
        df = pd.read_csv(INPUT_FILE, sep=r'\s+', header=None, names=['idx1', 'idx2', 'score'])
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Filter by distance
    df['dist'] = df['idx2'] - df['idx1']
    df_dist = df[df['dist'] >= MIN_BIN_DIST].copy()
    
    if df_dist.empty:
        print("No loops found after distance filtering.")
        header = ["#chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "color"]
        pd.DataFrame(columns=header).to_csv(OUTPUT_FILE, sep='\t', index=False)
        sys.exit(0)

    # Calculate the threshold by provided percentile 
    threshold = np.percentile(df_dist['score'], PERCENTILE)
    
    print(f"Calculated Threshold (Top {100-PERCENTILE:.1f}%): {threshold:.5f}")

    # Filter out loops by calculated threshold
    loops = df_dist[df_dist['score'] >= threshold].copy()
    print(f"Loops found: {len(loops)}")

    # Convert into coordinates
    loops['chr1'] = CHROM
    loops['x1'] = loops['idx1'] * RES
    loops['x2'] = loops['x1'] + RES 
    loops['chr2'] = CHROM
    loops['y1'] = loops['idx2'] * RES
    loops['y2'] = loops['y1'] + RES
    
    loops['name'] = '.' 
    loops['score_out'] = loops['score'].round(5)
    loops['color'] = '0,0,0' 

    bedpe_cols = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'name', 'score_out', 'color']
    header = ["#chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "color"]
    
    loops[bedpe_cols].to_csv(OUTPUT_FILE, sep='\t', index=False, header=header)
    print(f"Saved: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()