import pandas as pd
import numpy as np
import os
import sys
import argparse
import tempfile
import heapq


DEFAULT_CHUNK_SIZE = 2_000_000  # rows per chunk (~150 MB RAM per chunk)


def parse_args():
    parser = argparse.ArgumentParser(description="Process Juicer dump files (Observed & OE) for DeepLoop.")
    parser.add_argument("--obs",        required=True,                help="Path to 'dump observed' file")
    parser.add_argument("--oe",         required=True,                help="Path to 'dump oe' file")
    parser.add_argument("--chrom",      required=True,                help="Chromosome name")
    parser.add_argument("--res",        type=int, required=True,      help="Resolution in BP")
    parser.add_argument("--out",        required=True,                help="Output TSV filename")
    parser.add_argument("--anchor-dir", default="anchors",            help="Directory to save BED reference")
    parser.add_argument("--chunk-size", type=int,
                        default=DEFAULT_CHUNK_SIZE,
                        help=f"Rows per chunk for memory-efficient processing (Default: {DEFAULT_CHUNK_SIZE:,})")
    return parser.parse_args()


def detect_columns(filename):
    """Peek at first non-comment row to determine column layout."""
    with open(filename) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            n = len(line.split())
            if n >= 5:
                return [1, 3, 4]   # sparse full format: bin1 . bin2 . value
            elif n == 3:
                return [0, 1, 2]   # compact format: bin1 bin2 value
            else:
                raise ValueError(f"Unknown column layout ({n} cols) in {filename}")
    raise ValueError(f"Empty or header-only file: {filename}")


def sort_dump_to_tempfile(filename, col_indices, val_col_name, tmpdir, chunk_size):
    """
    Read dump file in chunks, normalise columns, sort each chunk,
    write sorted chunks to temp files. Returns list of temp file paths
    and max position seen (for BED generation).
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    chunk_files = []
    max_pos = 0

    reader = pd.read_csv(
        filename,
        sep=r'\s+',
        header=None,
        comment='#',
        on_bad_lines='skip',
        chunksize=chunk_size,
        dtype=str,
    )

    for i, chunk in enumerate(reader):
        chunk = chunk.iloc[:, col_indices].copy()
        chunk.columns = ['pos1', 'pos2', val_col_name]

        chunk['pos1'] = pd.to_numeric(chunk['pos1'], errors='coerce').fillna(0).astype(int)
        chunk['pos2'] = pd.to_numeric(chunk['pos2'], errors='coerce').fillna(0).astype(int)
        chunk[val_col_name] = pd.to_numeric(chunk[val_col_name], errors='coerce')
        chunk.dropna(inplace=True)

        chunk.sort_values(['pos1', 'pos2'], inplace=True)

        max_pos = max(max_pos, int(chunk['pos1'].max()), int(chunk['pos2'].max()))

        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.tsv', dir=tmpdir, delete=False
        )
        chunk.to_csv(tmp, sep='\t', index=False, header=True)
        tmp.close()
        chunk_files.append(tmp.name)
        print(f"  Sorted chunk {i+1} ({len(chunk):,} rows)")

    return chunk_files, max_pos


def iter_sorted_file(path, val_col):
    """Yield (pos1, pos2, value) from a sorted chunk file."""
    for chunk in pd.read_csv(path, sep='\t', chunksize=100_000,
                              dtype={'pos1': int, 'pos2': int}):
        for row in chunk.itertuples(index=False):
            yield int(row.pos1), int(row.pos2), float(getattr(row, val_col))


def merge_sorted_chunks(chunk_files, val_col):
    """k-way merge of sorted chunk files using heapq. Yields (pos1, pos2, value)."""
    iters = [iter_sorted_file(f, val_col) for f in chunk_files]
    heap  = []

    for idx, it in enumerate(iters):
        row = next(it, None)
        if row is not None:
            heapq.heappush(heap, (row[0], row[1], idx, row[2]))

    while heap:
        p1, p2, idx, val = heapq.heappop(heap)
        yield p1, p2, val
        nxt = next(iters[idx], None)
        if nxt is not None:
            heapq.heappush(heap, (nxt[0], nxt[1], idx, nxt[2]))


def stream_merge_obs_oe(obs_chunks, oe_chunks, res, out_path):
    """
    Stream-merge two sorted sequences (obs, oe) on (pos1, pos2).
    Write output directly to disk — no full DataFrame in memory.
    Returns (max_anchor, n_written).
    """
    obs_iter = merge_sorted_chunks(obs_chunks, 'observed')
    oe_iter  = merge_sorted_chunks(oe_chunks,  'oe')

    max_anchor = 0
    written    = 0

    obs_row = next(obs_iter, None)
    oe_row  = next(oe_iter,  None)

    with open(out_path, 'w') as fout:
        while obs_row is not None and oe_row is not None:
            op1, op2, obs_val = obs_row
            ep1, ep2, oe_val  = oe_row

            if (op1, op2) == (ep1, ep2):
                expected = obs_val / oe_val if oe_val != 0 else 0.0
                if not np.isfinite(expected):
                    expected = 0.0

                a1 = op1 // res
                a2 = op2 // res
                max_anchor = max(max_anchor, a1, a2)

                fout.write(f"{a1}\t{a2}\t{obs_val:.6g}\t{expected:.6g}\n")
                written += 1

                obs_row = next(obs_iter, None)
                oe_row  = next(oe_iter,  None)

            elif (op1, op2) < (ep1, ep2):
                obs_row = next(obs_iter, None)
            else:
                oe_row  = next(oe_iter,  None)

    return max_anchor, written


def write_bed_reference(anchor_dir, chrom, max_anchor, res, padding_bins=200):
    """Write BED anchor reference file for DeepLoop."""
    os.makedirs(anchor_dir, exist_ok=True)
    bed_path = os.path.join(anchor_dir, f"{chrom}.bed")
    with open(bed_path, 'w') as f:
        for i in range(int(max_anchor) + padding_bins):
            s = i * res
            e = s + res
            f.write(f"{chrom}\t{s}\t{e}\t{i}\n")
    print(f"BED file created: {bed_path} (padding: {padding_bins} bins)")
    return bed_path


def main():
    args = parse_args()

    print(f"Processing Juicer dump files for {args.chrom} at {args.res} bp...")
    print(f"Chunk size: {args.chunk_size:,} rows")

    try:
        obs_cols = detect_columns(args.obs)
        oe_cols  = detect_columns(args.oe)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)

    with tempfile.TemporaryDirectory() as tmpdir:
        print("Sorting OBS file in chunks...")
        try:
            obs_chunks, _ = sort_dump_to_tempfile(args.obs, obs_cols, 'observed', tmpdir, args.chunk_size)
        except (FileNotFoundError, ValueError, RuntimeError) as e:
            print(f"Error reading OBS: {e}")
            sys.exit(1)

        print("Sorting OE file in chunks...")
        try:
            oe_chunks, _ = sort_dump_to_tempfile(args.oe, oe_cols, 'oe', tmpdir, args.chunk_size)
        except (FileNotFoundError, ValueError, RuntimeError) as e:
            print(f"Error reading OE: {e}")
            sys.exit(1)

        print(f"Stream-merging -> {args.out} ...")
        try:
            max_anchor, n_written = stream_merge_obs_oe(
                obs_chunks, oe_chunks, args.res, args.out
            )
        except Exception as e:
            print(f"Error during merge: {e}")
            sys.exit(1)

    print(f"Interactions written: {n_written:,}")

    print(f"Generating BED reference in: {args.anchor_dir} ...")
    try:
        bed_path = write_bed_reference(args.anchor_dir, args.chrom, max_anchor, args.res)
    except OSError as e:
        print(f"Error writing BED: {e}")
        sys.exit(1)

    print("Data ready for DeepLoop.")


if __name__ == "__main__":
    main()