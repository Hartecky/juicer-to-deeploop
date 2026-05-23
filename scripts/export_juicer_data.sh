#!/bin/bash

set -e

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Export Observed and OE matrices from .hic file using Juicer Tools."
    echo ""
    echo "Required Arguments:"
    echo "  -j, --jar    <FILE>   Path to juicer_tools.jar"
    echo "  -i, --input  <FILE>   Path to .hic file"
    echo "  -c, --chrom  <STR>    Chromosome name (e.g., chr1)"
    echo "  -r, --res    <INT>    Resolution in bp (e.g., 10000)"
    echo "  -n, --norm   <STR>    Normalization type (e.g., KR)"
    echo "  -o, --out    <STR>    Output prefix (e.g., results/chr1_10k)"
    echo ""
}

JAR=""
INTER_HIC=""
CHROM=""
RES=""
NORM=""
OUT_PREFIX=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -j|--jar)   JAR="$2";       shift ;;
        -i|--input) INTER_HIC="$2"; shift ;;
        -c|--chrom) CHROM="$2";     shift ;;
        -r|--res)   RES="$2";       shift ;;
        -n|--norm)  NORM="$2";      shift ;;
        -o|--out)   OUT_PREFIX="$2"; shift ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "Unknown parameter: $1"; usage; exit 1 ;;
    esac
    shift
done

if [ -z "$JAR" ] || [ -z "$INTER_HIC" ] || [ -z "$CHROM" ] || [ -z "$RES" ] || [ -z "$NORM" ] || [ -z "$OUT_PREFIX" ]; then
    echo "Error: Missing required arguments."
    usage
    exit 1
fi

if [ ! -f "$INTER_HIC" ]; then
    echo "Error: .hic file not found: $INTER_HIC"
    exit 1
fi

if [ ! -f "$JAR" ]; then
    echo "Error: juicer_tools.jar not found: $JAR"
    exit 1
fi

echo "[1/2] Dumping Observed matrix..."
java -jar "$JAR" dump observed "$NORM" "$INTER_HIC" "$CHROM" "$CHROM" BP "$RES" "${OUT_PREFIX}_obs.txt"
echo "Created: ${OUT_PREFIX}_obs.txt"

echo "[2/2] Dumping OE matrix..."
java -jar "$JAR" dump oe "$NORM" "$INTER_HIC" "$CHROM" "$CHROM" BP "$RES" "${OUT_PREFIX}_oe.txt"
echo "Created: ${OUT_PREFIX}_oe.txt"