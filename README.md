# juicer-to-deeploop
Bridge between juicer's .hic output file to DeepLoop algorithm for peak detection from Hi-C data

**juicer-to-deeploop** is a streamlined pipeline that automates the detection of chromatin loops from Hi-C data. It acts as a bridge between **Juicer Tools** (for data extraction) and **DeepLoop** (for deep learning based denoising), producing loop calls in `.bedpe` format compatible with Juicebox for visualization with actual .hic files from juicer pipeline.

# Pipeline Logic
1. Juicer Dump: Extracts "Observed" (Raw) and "Observed/Expected" (OE) KR-normalized matrices from the .hic file.
2. Pre-processing: Merges matrices, calculates expected values (required for DeepLoop), and handles padding to ensure DeepLoop doesn't crash on chromosome edges.
3. DeepLoop: Runs the neural network to denoise the map and assign probabilities to pixels.
4. Filtering: Converts the probability map to a list of loops, filtering out the diagonal (short-range artifacts) and applying a percentile threshold.

## Prerequisites
Before running the pipeline, ensure you have the following installed:

1. Juicer - developed by Aidenlab

   https://github.com/aidenlab/juicer/tree/main?tab=readme-ov-file
   
   Neva C. Durand, Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments." Cell Systems 3(1), 2016
   
2. DeepLoop - developed by JinLab

   https://github.com/JinLabBioinfo/DeepLoop

   Zhang, S., Plummer, D., Lu, L. et al. DeepLoop robustly maps chromatin interactions from sparse allele-resolved or single-cell Hi-C data at kilobase resolution. Nat Genet 54, 1013–1025 (2022)

## Repository Structure
```text
├── README.md                     
├── config/
│   └── config.yaml               
├── envs/
│   ├── environment.yaml          
│   └── requirements.txt         
├── juicer_deeploop_pipeline.sh   
└── scripts/
    ├── deeploop_to_bedpe.py      
    ├── export_juicer_data.sh     
    └── process_juicer_output.py
```

## Installation
### Clone the repository
```bash
git clone https://github.com/Hartecky/juicer-to-deeploop
cd juicer-to-deeploop
```

### Setup environments
```bash
# Create the Conda environment from the YAML file
conda env create -f envs/environment.yaml

# Activate the environment
conda activate juicer-to-deeploop-env

# Install python dependencies
pip install -r envs/requirements.txt
```
### Setup paths in workflow scripts

**UPDATE THESE PATHS IN juicer_to_deeploop_pipeline.sh SCRIPT**

Note: The pipeline assumes that pre-trained DeepLoop models are in $DL_DIR/DeepLoop_models/CPGZ_trained/ downloaded from JinLab repository

```bash
JUICER_JAR="/absolute/path/to/juicer_tools.jar"
DL_DIR="/absolute/path/to/DeepLoop" 
```
### Optional tuning parameters
Inside juicer_deeploop_pipeline.sh, you can also adjust the loop filtering thresholds:

  -  FIXED_MIN_DIST=5: Minimum distance in bins from the diagonal (removes local noise).
  -  FIXED_PERCENTILE=98.0: Keeps only the top 2% of strongest signals detected by DeepLoop.

## Usage

Run the pipeline using the main bash script

```bash
./juicer_deeploop_pipeline.sh -i <hic_file> -c <chrom> -r <resolution> -o <output_dir>

    -i, --input: Path to the .hic file (e.g., inter_30.hic)
    -c, --chrom: Chromosome name (must match .hic internals, usually chr1, chr2)
    -r, --res: Resolution in base pairs (e.g., 2000, 5000, 10000)
    -o, --out: Output directory for results
```

Example
```bash
./juicer_deeploop_pipeline.sh \
  --input data/GM12878.hic \
  --chrom chr1 \
  --res 5000 \
  --out results/chr1_5k
```

## Outputs
The pipeline creates the following directory structure in your output folder:

  -  final_bedpe/: Main Result. Contains .bedpe files with detected loops (loadable in Juicebox).
  -  deeploop_out/: Raw output from DeepLoop (denoised probability matrices).
  -  deeploop_in/: Input text files prepared for DeepLoop.
  -  raw_dumps/: Raw Observed and OE matrices dumped from Juicer.
  -  anchors/: Genomic coordinates (BED files) used for mapping bins.
