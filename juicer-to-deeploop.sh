#!/bin/bash
set -e

# Help function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Multiscale Pipeline for Chromatin Loop Detection (DeepLoop + DBSCAN)."
    echo ""
    echo "Required Arguments:"
    echo "  -i, --input      <FILE>    Path to input .hic file"
    echo "  -c, --chrom      <STR>     Chromosome name (e.g., chr1)"
    echo "  -r, --res        <LIST>    Comma-separated resolutions (e.g., 2000,5000,10000)"
    echo "  -n, --norm       <STR>     Normalization type (one of KR, SCALE, VC, VC_SQRT, GW_SCALE, INTER_SCALE)"
    echo "  -o, --out        <DIR>     Output directory path"
    echo "  -t, --tolerance  <INT>     Merge tolerance in bp (Default: 20000)"
    echo "  -k, --chunk-size <INT>     Rows per chunk for process_juicer_output.py (Default: 2000000)"
    echo ""
    echo "Example:"
    echo "  $0 --input data/inter.hic --chrom chr1 --res 2000,5000,10000 --norm KR --out results_chr1"
    echo ""
}

# Parse Arguments
HIC_FILE=""
CHROM=""
RES_STRING=""
NORM_STRING=""
OUT_DIR=""
MERGE_TOLERANCE=20000
CHUNK_SIZE=2000000

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input)      HIC_FILE="$2";       shift ;;
        -c|--chrom)      CHROM="$2";          shift ;;
        -r|--res)        RES_STRING="$2";     shift ;;
        -n|--norm)       NORM_STRING="$2";    shift ;;
        -o|--out)        OUT_DIR="$2";        shift ;;
        -t|--tolerance)  MERGE_TOLERANCE="$2"; shift ;;
        -k|--chunk-size) CHUNK_SIZE="$2";     shift ;;
        -h|--help)       usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
done

if [ -z "$HIC_FILE" ] || [ -z "$CHROM" ] || [ -z "$RES_STRING" ] || [ -z "$NORM_STRING" ] || [ -z "$OUT_DIR" ]; then
    echo "Error: Missing required arguments."
    usage
    exit 1
fi

# Load config
CONFIG_FILE="$(dirname "$0")/config/config.yaml"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: config.yaml not found at $CONFIG_FILE"
    exit 1
fi

parse_yaml() {
    local key="$1"
    grep "^${key}:" "$CONFIG_FILE" \
        | sed 's/^[^:]*: *//' \
        | tr -d '"' \
        | sed "s|\${HOME}|$HOME|g"
}

JUICER_JAR="$(parse_yaml juicer_jar)"
DL_DIR="$(parse_yaml deeploop_dir)"
MODEL_H5="$DL_DIR/$(parse_yaml model_h5)"
MODEL_JSON="$DL_DIR/$(parse_yaml model_json)"

# Validate paths from config
if [ ! -f "$JUICER_JAR" ]; then
    echo "Error: juicer_tools.jar not found at: $JUICER_JAR"
    echo "Update juicer_jar in config/config.yaml"
    exit 1
fi

if [ ! -d "$DL_DIR" ]; then
    echo "Error: DeepLoop directory not found at: $DL_DIR"
    echo "Update deeploop_dir in config/config.yaml"
    exit 1
fi

SCRIPT_PROCESS="scripts/process_juicer_output.py"
SCRIPT_DBSCAN="scripts/deeploop_to_bedpe.py"
SCRIPT_MERGE="scripts/merge_loops.py"
SCRIPT_DUMP="scripts/export_juicer_data.sh"

# Create directories
mkdir -p "$OUT_DIR"/{raw_dumps,anchors,deeploop_in,deeploop_out,final_bedpe}

# Setup logging
LOG_DIR="$OUT_DIR/logs"
mkdir -p "$LOG_DIR"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="$LOG_DIR/pipeline_${CHROM}_${TIMESTAMP}.log"

# Tee all output to log file
exec > >(tee -a "$LOG_FILE") 2>&1

PIPELINE_START=$(date +%s)
echo "========================================"
echo "START MULTISCALE PIPELINE"
echo "Date:        $(date '+%Y-%m-%d %H:%M:%S')"
echo "Chromosome:  $CHROM"
echo "Resolutions: $RES_STRING"
echo "Norm:        $NORM_STRING"
echo "Input:       $HIC_FILE"
echo "Output dir:  $OUT_DIR"
echo "Tolerance:   $MERGE_TOLERANCE bp"
echo "Chunk size:  $CHUNK_SIZE"
echo "Log file:    $LOG_FILE"
echo "========================================"

# Split resolutions string into array (comma delimiter)
IFS=',' read -ra RES_ARRAY <<< "$RES_STRING"

# Arrays to store paths and loop counts for the final merge and summary
FILES_TO_MERGE=()
RES_TO_MERGE=()
declare -A LOOPS_PER_RES

# Iterate over resolutions
for RES in "${RES_ARRAY[@]}"; do
    echo ""
    echo "----------------------------------------"
    echo "Processing resolution: $RES bp"
    echo "----------------------------------------"
    RES_START=$(date +%s)

    PREFIX="$OUT_DIR/raw_dumps/${CHROM}_${RES}"
    DL_INPUT="$OUT_DIR/deeploop_in/${CHROM}_${RES}_input.txt"
    DL_OUTPUT="$OUT_DIR/deeploop_out/${CHROM}_${RES}.denoised.anchor.to.anchor"
    BEDPE_OUTPUT="$OUT_DIR/final_bedpe/${CHROM}_${RES}_loops.bedpe"

    # CONFIGURATION OF SENSITIVITY
    if (( RES == 25000 )); then
        CURRENT_THRESHOLD=0.70
        CURRENT_MIN_DIST=5
        CURRENT_EPS=2.5
        CURRENT_MIN_SAMPLES=2
    elif (( RES == 10000 )); then
        CURRENT_THRESHOLD=0.70
        CURRENT_MIN_DIST=5
        CURRENT_EPS=2.5
        CURRENT_MIN_SAMPLES=2
    elif (( RES == 5000 )); then
        CURRENT_THRESHOLD=0.80
        CURRENT_MIN_DIST=4
        CURRENT_EPS=3.0
        CURRENT_MIN_SAMPLES=2
    elif (( RES == 2000 )); then
        CURRENT_THRESHOLD=0.90
        CURRENT_MIN_DIST=3
        CURRENT_EPS=3.5
        CURRENT_MIN_SAMPLES=2
    elif (( RES == 1000 )); then
        CURRENT_THRESHOLD=0.90
        CURRENT_MIN_DIST=3
        CURRENT_EPS=3.5
        CURRENT_MIN_SAMPLES=2
    elif (( RES == 500 )); then
        CURRENT_THRESHOLD=0.90
        CURRENT_MIN_DIST=3
        CURRENT_EPS=3.5
        CURRENT_MIN_SAMPLES=2
    else
        echo "Unsupported resolution: RES=$RES" >&2
        exit 1
    fi

    # [1/5] Juicer dump
    echo "[1/5] Dumping data..."
    if [ ! -f "${PREFIX}_obs.txt" ]; then
        bash "$SCRIPT_DUMP" \
            --jar   "$JUICER_JAR" \
            --input "$HIC_FILE" \
            --chrom "$CHROM" \
            --res   "$RES" \
            --norm  "$NORM_STRING" \
            --out   "$PREFIX"
    else
        echo "Dump exists, skipping."
    fi

    # [2/5] Prepare input for DeepLoop
    echo "[2/5] Processing input..."
    python "$SCRIPT_PROCESS" \
        --obs        "${PREFIX}_obs.txt" \
        --oe         "${PREFIX}_oe.txt" \
        --chrom      "$CHROM" \
        --res        "$RES" \
        --out        "$DL_INPUT" \
        --anchor-dir "$OUT_DIR/anchors" \
        --chunk-size "$CHUNK_SIZE"

    # [3/5] Run DeepLoop model
    echo "[3/5] Running DeepLoop..."
    export CUDA_VISIBLE_DEVICES=""

    RES_OUT_DIR="$OUT_DIR/deeploop_out/$RES"
    mkdir -p "$RES_OUT_DIR"

    if [ ! -f "$DL_OUTPUT" ]; then
        python "$DL_DIR/prediction/predict_chromosome.py" \
            --full_matrix_dir "$OUT_DIR/deeploop_in" \
            --input_name      "$(basename "$DL_INPUT")" \
            --out_dir         "$RES_OUT_DIR" \
            --anchor_dir      "$OUT_DIR/anchors" \
            --h5_file         "$MODEL_H5" \
            --json_file       "$MODEL_JSON" \
            --chromosome      "$CHROM" \
            --val_cols obs exp \
            --small_matrix_size 128 \
            --step_size 64 \
            --dummy 5 > /dev/null

        mv "$RES_OUT_DIR/${CHROM}.denoised.anchor.to.anchor" "$DL_OUTPUT"
    else
        echo "DeepLoop output exists, skipping prediction."
    fi

    if [ ! -f "$DL_OUTPUT" ]; then
        echo "ERROR: DeepLoop failed for $RES." >&2
        exit 1
    fi

    # [4/5] DBSCAN clustering -> BEDPE
    echo "[4/5] Clustering (DBSCAN) -> BEDPE..."
    python "$SCRIPT_DBSCAN" \
        --input      "$DL_OUTPUT" \
        --out        "$BEDPE_OUTPUT" \
        --chrom      "$CHROM" \
        --res        "$RES" \
        --min-dist   "$CURRENT_MIN_DIST" \
        --threshold  "$CURRENT_THRESHOLD" \
        --eps        "$CURRENT_EPS" \
        --min-samples "$CURRENT_MIN_SAMPLES"

    # Count loops and store for summary
    if [ -f "$BEDPE_OUTPUT" ]; then
        N_LOOPS=$(tail -n +2 "$BEDPE_OUTPUT" | wc -l)
        LOOPS_PER_RES[$RES]=$N_LOOPS
        echo "Loops detected at ${RES} bp: $N_LOOPS"
        FILES_TO_MERGE+=("$BEDPE_OUTPUT")
        RES_TO_MERGE+=("$RES")
    fi

    RES_END=$(date +%s)
    echo "Resolution $RES done in $(( RES_END - RES_START ))s"

done

# [5/5] Merge multiresolution results
echo ""
echo "----------------------------------------"
echo "[5/5] Merging multiresolution results..."
echo "----------------------------------------"

FINAL_MERGED="$OUT_DIR/final_bedpe/${CHROM}_merged_multires.bedpe"

python "$SCRIPT_MERGE" \
    --files       "${FILES_TO_MERGE[@]}" \
    --resolutions "${RES_TO_MERGE[@]}" \
    --out         "$FINAL_MERGED" \
    --tolerance   "$MERGE_TOLERANCE"

N_MERGED=$(tail -n +2 "$FINAL_MERGED" | wc -l)

PIPELINE_END=$(date +%s)
PIPELINE_ELAPSED=$(( PIPELINE_END - PIPELINE_START ))

# Write final summary
SUMMARY_FILE="$LOG_DIR/final_summary_${CHROM}_${TIMESTAMP}.txt"
{
    echo "========================================"
    echo "PIPELINE SUMMARY"
    echo "========================================"
    echo "Date:           $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Chromosome:     $CHROM"
    echo "Resolutions:    $RES_STRING"
    echo "Normalization:  $NORM_STRING"
    echo "Input:          $HIC_FILE"
    echo "Output dir:     $OUT_DIR"
    echo "Tolerance:      $MERGE_TOLERANCE bp"
    echo "Chunk size:     $CHUNK_SIZE"
    echo ""
    echo "--- Loops per resolution ---"
    for RES in "${RES_ARRAY[@]}"; do
        printf "  %6d bp:  %d loops\n" "$RES" "${LOOPS_PER_RES[$RES]:-0}"
    done
    echo ""
    echo "--- Final merged result ---"
    echo "  Loops after merge:  $N_MERGED"
    echo "  Output file:        $FINAL_MERGED"
    echo ""
    echo "Total runtime:  ${PIPELINE_ELAPSED}s ($(( PIPELINE_ELAPSED / 60 ))m $(( PIPELINE_ELAPSED % 60 ))s)"
    echo "Log file:       $LOG_FILE"
    echo "========================================"
} | tee "$SUMMARY_FILE"

echo ""
echo "DONE. Summary saved to: $SUMMARY_FILE"