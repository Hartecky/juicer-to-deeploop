import pandas as pd
import numpy as np
import sys
import os
import argparse
from sklearn.cluster import DBSCAN

BEDPE_HEADER = ["#chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "color"]


def parse_args():
    parser = argparse.ArgumentParser(description="Convert DeepLoop output to BEDPE using DBSCAN clustering.")
    parser.add_argument("--input", required=True, help="Path to DeepLoop output file")
    parser.add_argument("--out", required=True, help="Output BEDPE filename")
    parser.add_argument("--chrom", required=True, help="Chromosome name")
    parser.add_argument("--res", type=int, required=True, help="Resolution in BP")

    # DBSCAN parameters
    parser.add_argument("--min-dist", type=int, default=5, help="Min distance in bins (Default: 5)")
    parser.add_argument("--threshold", type=float, default=0.97, help="Score threshold for DBSCAN input (Default: 0.97)")
    parser.add_argument("--eps", type=float, default=2.5, help="DBSCAN eps (in bins). Default: 2.5")
    parser.add_argument("--min-samples", type=int, default=3, help="DBSCAN min_samples. Default: 3")

    return parser.parse_args()


def write_empty_bedpe(path):
    """Write an empty BEDPE file with correct header and exit."""
    pd.DataFrame(columns=BEDPE_HEADER).to_csv(path, sep='\t', index=False)
    sys.exit(0)


def main():

    args = parse_args()

    print(f"--- DeepLoop DBSCAN: {args.chrom} @ {args.res}bp ---")

    if not os.path.exists(args.input):
        print(f"Error: Missing file {args.input}")
        sys.exit(1)

    try:
        df = pd.read_csv(args.input, sep=r'\s+', header=None, names=['idx1', 'idx2', 'score'])
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    # 1. Pre-filtering
    df['dist'] = df['idx2'] - df['idx1']
    # Remove diagonal (self-interactions)
    df = df[df['dist'] >= args.min_dist]
    # Remove low-confidence interactions below score threshold
    df = df[df['score'] >= args.threshold].copy()

    if df.empty:
        print("No loops found after filtering.")
        write_empty_bedpe(args.out)

    # 2. DBSCAN clustering
    print(f"Running DBSCAN on {len(df)} points...")
    coords = df[['idx1', 'idx2']].values

    # eps defines the max distance (in bins) between points in the same cluster
    db = DBSCAN(eps=args.eps, min_samples=args.min_samples).fit(coords)

    df['cluster'] = db.labels_

    # Discard noise points (label == -1)
    clustered = df[df['cluster'] != -1]

    unique_clusters = clustered['cluster'].nunique()
    print(f"Found {unique_clusters} loops (clusters).")

    if unique_clusters == 0:
        print("No clusters found after DBSCAN.")
        write_empty_bedpe(args.out)

    # 3. Compute cluster centroids
    loops_data = []
    for cid, group in clustered.groupby('cluster'):
        # Center of mass
        c_idx1 = int(group['idx1'].mean())
        c_idx2 = int(group['idx2'].mean())
        # Peak score within the cluster
        c_score = group['score'].max()

        loops_data.append([c_idx1, c_idx2, c_score])

    final_df = pd.DataFrame(loops_data, columns=['idx1', 'idx2', 'score'])

    # 4. Convert bin indices to genomic coordinates
    final_df['chr1'] = args.chrom
    final_df['x1'] = final_df['idx1'] * args.res
    final_df['x2'] = final_df['x1'] + args.res
    final_df['chr2'] = args.chrom
    final_df['y1'] = final_df['idx2'] * args.res
    final_df['y2'] = final_df['y1'] + args.res

    final_df['name'] = '.'
    final_df['score_out'] = final_df['score'].round(5)
    final_df['color'] = '0,0,0'

    cols = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'name', 'score_out', 'color']

    final_df[cols].to_csv(args.out, sep='\t', index=False, header=BEDPE_HEADER)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()