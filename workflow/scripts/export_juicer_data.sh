#!/bin/bash

set -e

# Help function with instructions about usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Export KR normalized Observed and OE matrices from .hic file using Juicer Tools."
    echo ""
    echo "Required Arguments:"
    echo "  -i, --input <FILE>      Path to .hic file"
    echo "  -c, --chrom <STR>       Chromosome name (e.g., chr1)"
    echo "  -r, --res <INT>         Resolution in bp (e.g., 10000)"
    echo "  -o, --out-prefix <STR>  Output filename prefix (e.g., output_chr1)"
    echo ""
    echo "Optional Arguments:"
    echo "  -h, --help              Show this help message and exit"
    echo ""
    echo "Example:"
    echo "  $0 --input ../aligned/inter_30.hic --chrom chr1 --res 10000 --out-prefix results/chr1_data"
    echo ""
}

# Arguments parsing - initialize variables as empty vectors which will be filled from user
INTER_HIC=""
CHROM=""
RES=""
OUT_PREFIX=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INTER_HIC="$2"; shift ;;
        -c|--chrom) CHROM="$2"; shift ;;
        -r|--res) RES="$2"; shift ;;
        -o|--out-prefix) OUT_PREFIX="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
done

# Arguments validation
if [ -z "$INTER_HIC" ] || [ -z "$CHROM" ] || [ -z "$RES" ] || [ -z "$OUT_PREFIX" ]; then
    echo "Error: Missing required arguments."
    echo "---------------------------------"
    usage
    exit 1
fi


JUICER_JAR="/mnt/storage_6/project_data/pl0457-01/pekowska_lab_software/juicer/scripts/common/juicer_tools.jar"

if [ ! -f "$INTER_HIC" ]; then
    echo "Error: Hi-C file does not exist: $INTER_HIC"
    exit 1
fi

if [ ! -f "$JUICER_JAR" ]; then
    echo "Error: juicer_tools.jar not found at: $JUICER_JAR"
    exit 1
fi

# Dump Observed data (Raw/Observed Matrix)
echo "[1/2] Dumping Observed matrix..."
java -jar "$JUICER_JAR" dump observed KR "$INTER_HIC" "$CHROM" "$CHROM" BP "$RES" "${OUT_PREFIX}_obs.txt"
echo "Created: ${OUT_PREFIX}_obs.txt"

# Dump OE data (Observed/Expected)
echo "[2/2] Dumping OE matrix..."
java -jar "$JUICER_JAR" dump oe KR "$INTER_HIC" "$CHROM" "$CHROM" BP "$RES" "${OUT_PREFIX}_oe.txt"
echo "Created: ${OUT_PREFIX}_oe.txt"