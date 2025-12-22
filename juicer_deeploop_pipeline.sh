#!/bin/bash

# Help function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Multiscale Pipeline for Chromatin Loop Detection (DeepLoop + DBSCAN)."
    echo ""
    echo "Required Arguments:"
    echo "  -i, --input <FILE>    Path to input .hic file"
    echo "  -c, --chrom <STR>     Chromosome name (e.g., chr1)"
    echo "  -r, --res <LIST>      Comma-separated resolutions (e.g., 2000,5000,10000)"
    echo "  -o, --out <DIR>       Output directory path"
    echo ""
    echo "Example:"
    echo "  $0 --input data/inter.hic --chrom chr1 --res 2000,5000,10000 --out results_chr1"
    echo ""
}

# Parse Arguments
HIC_FILE=""
CHROM=""
RES_STRING="" # String np "2000,5000"
OUT_DIR=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) HIC_FILE="$2"; shift ;;
        -c|--chrom) CHROM="$2"; shift ;;
        -r|--res) RES_STRING="$2"; shift ;;
        -o|--out) OUT_DIR="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
done

if [ -z "$HIC_FILE" ] || [ -z "$CHROM" ] || [ -z "$RES_STRING" ] || [ -z "$OUT_DIR" ]; then
    echo "Error: Missing arguments."
    usage
    exit 1
fi

# Configs (Update Paths!)
JUICER_JAR="/mnt/storage_6/project_data/pl0457-01/pekowska_lab_software/juicer/scripts/common/juicer_tools.jar"
DL_DIR="/mnt/storage_3/home/b.hofman/pl0457-01/project_data/pekowska_lab_software/DeepLoop"
MODEL_H5="$DL_DIR/DeepLoop_models/CPGZ_trained/LoopDenoise.h5"
MODEL_JSON="$DL_DIR/DeepLoop_models/CPGZ_trained/LoopDenoise.json"

SCRIPT_PROCESS="scripts/process_juicer_output.py"
SCRIPT_DBSCAN="scripts/deeploop_to_bedpe.py" # NOWY SKRYPT
SCRIPT_MERGE="scripts/merge_loops.py" # NOWY SKRYPT

# Create directories
mkdir -p "$OUT_DIR/raw_dumps"
mkdir -p "$OUT_DIR/anchors"
mkdir -p "$OUT_DIR/deeploop_in"
mkdir -p "$OUT_DIR/deeploop_out"
mkdir -p "$OUT_DIR/final_bedpe"

echo "START MULTISCALE PIPELINE: $CHROM"
echo "Resolutions: $RES_STRING"

# Split resolutions string into array (comma delimiter)
IFS=',' read -ra RES_ARRAY <<< "$RES_STRING"

# Arrays to store paths for the final merge step
FILES_TO_MERGE=()
RES_TO_MERGE=()

# Iterate over three resolutions
for RES in "${RES_ARRAY[@]}"; do
    echo ""
    echo ">>> Processing resolution: $RES bp"
    
    PREFIX="$OUT_DIR/raw_dumps/${CHROM}_${RES}"
    DL_INPUT="$OUT_DIR/deeploop_in/${CHROM}_${RES}_input.txt"
    DL_OUTPUT="$OUT_DIR/deeploop_out/${CHROM}_${RES}.denoised.anchor.to.anchor" # Dodalem RES do nazwy
    BEDPE_OUTPUT="$OUT_DIR/final_bedpe/${CHROM}_${RES}_loops.bedpe"

    # 1. Juicer Dump
    echo "[1/4] Dumping data..."
    if [ ! -f "${PREFIX}_obs.txt" ]; then
        java -jar "$JUICER_JAR" dump observed KR "$HIC_FILE" "$CHROM" "$CHROM" BP "$RES" "${PREFIX}_obs.txt"
        java -jar "$JUICER_JAR" dump oe KR "$HIC_FILE" "$CHROM" "$CHROM" BP "$RES" "${PREFIX}_oe.txt"
    else
        echo "Dump exists, skipping."
    fi

    # 2. Process for DeepLoop
    echo "[2/4] Processing input..."
    python "$SCRIPT_PROCESS" \
        --obs "${PREFIX}_obs.txt" \
        --oe "${PREFIX}_oe.txt" \
        --chrom "$CHROM" \
        --res "$RES" \
        --out "$DL_INPUT" \
        --anchor-dir "$OUT_DIR/anchors"

    # 3. DeepLoop (Warning: DeepLoop output filename is tricky to control)
    echo "[3/4] Running DeepLoop..."
    export CUDA_VISIBLE_DEVICES=""
    
    # DeepLoop zawsze zapisuje jako {CHROM}.denoised.anchor.to.anchor w out_dir
    # Żeby nie nadpisywać plików przy różnych rozdzielczościach, zrobimy trick z podkatalogiem
    RES_OUT_DIR="$OUT_DIR/deeploop_out/$RES"
    mkdir -p "$RES_OUT_DIR"

    python "$DL_DIR/prediction/predict_chromosome.py" \
      --full_matrix_dir "$OUT_DIR/deeploop_in" \
      --input_name "$(basename $DL_INPUT)" \
      --out_dir "$RES_OUT_DIR" \
      --anchor_dir "$OUT_DIR/anchors" \
      --h5_file "$MODEL_H5" \
      --json_file "$MODEL_JSON" \
      --chromosome "$CHROM" \
      --val_cols obs exp \
      --small_matrix_size 128 \
      --step_size 64 \
      --dummy 5 > /dev/null

    # Przenosimy wynik do głównego folderu ze zmienioną nazwą
    mv "$RES_OUT_DIR/${CHROM}.denoised.anchor.to.anchor" "$DL_OUTPUT"
    
    if [ ! -f "$DL_OUTPUT" ]; then
        echo "ERROR: DeepLoop failed for $RES."
        exit 1
    fi

    # 4. Convert to BEDPE using DBSCAN
    echo "[4/4] Clustering (DBSCAN) -> BEDPE..."
    
    # Parametry zależne od rozdzielczości (opcjonalnie można sterować)
    # Dla 2k/5k/10k stałe eps=2.5 bina zazwyczaj działa dobrze
    python "$SCRIPT_DBSCAN" \
        --input "$DL_OUTPUT" \
        --out "$BEDPE_OUTPUT" \
        --chrom "$CHROM" \
        --res "$RES" \
        --min-dist 5 \
        --threshold 0.97 \
        --eps 2.5 \
        --min-samples 3

    # Add to list for merging
    if [ -f "$BEDPE_OUTPUT" ]; then
        FILES_TO_MERGE+=("$BEDPE_OUTPUT")
        RES_TO_MERGE+=("$RES")
    fi

done

echo "MERGING MULTISCALE RESULTS"

FINAL_MERGED="$OUT_DIR/final_bedpe/${CHROM}_merged_multires.bedpe"

# Przekazujemy tablice jako listę argumentów
python "$SCRIPT_MERGE" \
    --files "${FILES_TO_MERGE[@]}" \
    --resolutions "${RES_TO_MERGE[@]}" \
    --out "$FINAL_MERGED" \
    --tolerance 20000

echo "DONE."
echo "Final file: $FINAL_MERGED"