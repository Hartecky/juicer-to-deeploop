#!/bin/bash

# Help function with instructions about usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Pipeline for Chromatin Loop Detection using DeepLoop algorithm."
    echo ""
    echo "Required Arguments:"
    echo "  -i, --input <FILE>    Path to input .hic file"
    echo "  -c, --chrom <STR>     Chromosome name (e.g., chr1)"
    echo "  -r, --res <INT>       Resolution in bp (e.g., 2000)"
    echo "  -o, --out <DIR>       Output directory path"
    echo ""
    echo "Optional Arguments:"
    echo "  -h, --help            Show this help message and exit"
    echo ""
    echo "Example:"
    echo "  $0 --input data/inter.hic --chrom chr1 --res 2000 --out results_chr1"
    echo ""
}

# Arguments parsing - initialize variables as empty vectors which will be filled from user
HIC_FILE=""
CHROM=""
RES=""
OUT_DIR=""

# While loop that captures the arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) HIC_FILE="$2"; shift ;;
        -c|--chrom) CHROM="$2"; shift ;;
        -r|--res) RES="$2"; shift ;;
        -o|--out) OUT_DIR="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
done

# Validation
# Sprawdzamy, czy wszystkie wymagane zmienne zostały ustawione
if [ -z "$HIC_FILE" ] || [ -z "$CHROM" ] || [ -z "$RES" ] || [ -z "$OUT_DIR" ]; then
    echo "Error: Missing required arguments. \n"
    # echo ""
    usage
    exit 1
fi

# Configs & paths to softwares
JUICER_JAR="/mnt/storage_6/project_data/pl0457-01/pekowska_lab_software/juicer/scripts/common/juicer_tools.jar"
DL_DIR="/mnt/storage_3/home/b.hofman/pl0457-01/project_data/pekowska_lab_software/DeepLoop"
MODEL_H5="$DL_DIR/DeepLoop_models/CPGZ_trained/LoopDenoise.h5"
MODEL_JSON="$DL_DIR/DeepLoop_models/CPGZ_trained/LoopDenoise.json"

# Paths to python scripts
SCRIPT_PROCESS_JUICER="/mnt/storage_3/home/b.hofman/pl0457-01/scratch/SD071125_XXXXXXH/deep_loop_pipeline/scripts/process_juicer_output.py"
SCRIPT_TO_BEDPE="/mnt/storage_3/home/b.hofman/pl0457-01/scratch/SD071125_XXXXXXH/deep_loop_pipeline/scripts/deeploop_to_bedpe.py"

# Filtration parameters for final detected loops
FIXED_MIN_DIST=5        # minimum distance in bins (e.g. 5 bins from diagonal)
FIXED_PERCENTILE=98.0    # Percentile (e.g. 98.0 means that we take top 2% of the strongest signals)

# Create a directory structure 
mkdir -p "$OUT_DIR/raw_dumps"
mkdir -p "$OUT_DIR/anchors"
mkdir -p "$OUT_DIR/deeploop_in"
mkdir -p "$OUT_DIR/deeploop_out"
mkdir -p "$OUT_DIR/final_bedpe"

PREFIX="$OUT_DIR/raw_dumps/${CHROM}_${RES}"
DL_INPUT="$OUT_DIR/deeploop_in/${CHROM}_${RES}_input.txt"

echo "=========================================================="
echo "START PIPELINE: $CHROM @ $RES bp"
echo "Output Directory: $OUT_DIR"
echo "Params: MinDist=${FIXED_MIN_DIST}, Percentile=${FIXED_PERCENTILE}%"
echo "=========================================================="


# Juicer Dump the observed and observed/expected values after KR normalization
echo "[1/4] Dumping data from juicer..."
java -jar "$JUICER_JAR" dump observed KR "$HIC_FILE" "$CHROM" "$CHROM" BP "$RES" "${PREFIX}_obs.txt"
java -jar "$JUICER_JAR" dump oe KR "$HIC_FILE" "$CHROM" "$CHROM" BP "$RES" "${PREFIX}_oe.txt"


# Processing juicer dump output
echo "[2/4] Processing juicer output for DeepLoop..."
python "$SCRIPT_PROCESS_JUICER" \
    --obs "${PREFIX}_obs.txt" \
    --oe "${PREFIX}_oe.txt" \
    --chrom "$CHROM" \
    --res "$RES" \
    --out "$DL_INPUT" \
    --anchor-dir "$OUT_DIR/anchors"

# Run DeepLoop
echo "[3/4] Running DeepLoop loop calling ..."
export CUDA_VISIBLE_DEVICES=""
python "$DL_DIR/prediction/predict_chromosome.py" \
  --full_matrix_dir "$OUT_DIR/deeploop_in" \
  --input_name "$(basename $DL_INPUT)" \
  --out_dir "$OUT_DIR/deeploop_out" \
  --anchor_dir "$OUT_DIR/anchors" \
  --h5_file "$MODEL_H5" \
  --json_file "$MODEL_JSON" \
  --chromosome "$CHROM" \
  --val_cols obs exp \
  --small_matrix_size 128 \
  --step_size 64 \
  --dummy 5

DL_OUTPUT="$OUT_DIR/deeploop_out/${CHROM}.denoised.anchor.to.anchor"

if [ ! -f "$DL_OUTPUT" ]; then
    echo "ERROR: DeepLoop output missing: $DL_OUTPUT"
    exit 1
fi


# Convert Deeploop output signals to BEDPE for compatibility with juicebox
echo "[4/4] Converting DeepLoop signals into BEDPE..."
FINAL_BEDPE="$OUT_DIR/final_bedpe/${CHROM}_${RES}_loops.bedpe"

python "$SCRIPT_TO_BEDPE" \
    --input "$DL_OUTPUT" \
    --out "$FINAL_BEDPE" \
    --chrom "$CHROM" \
    --res "$RES" \
    --min-dist "$FIXED_MIN_DIST" \
    --percentile "$FIXED_PERCENTILE"