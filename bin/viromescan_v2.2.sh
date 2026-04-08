#!/bin/bash
# viromescan v2.2
# Pipeline for virome profiling from shotgun metagenomic data
# Uses Bowtie2 for host/microbiome filtering and Salmon EM for quantification
# Author: Hamidreza Chehrehnejad — MSc Molecular Biology and Genetics
# Thesis: Development and validation of a new pipeline for virome profiling

set -e

if [ $# -eq 0 ]; then
echo -e "\n$0 -1 <R1> [-2 <R2>] --host <HOST> -o <OUTPUT> [OPTIONS]

Required:
  -1/--input1       R1 (or single-end) FASTQ [.gz accepted]
  -2/--input2       R2 FASTQ (paired-end only)
  --host            cattle | pig | grapevine | human
  -o/--output       output directory

Optional:
  -p/--threads      threads (default: 8)
  --minreads        min Salmon NumReads for detection (default: 50)
  --mintpm          min TPM for detection (default: 1.0)
  --salmon-index    path to Salmon index (overrides default)
"
exit 1
fi

#  Argument parsing 
INPUT1=''; INPUT2=''; HOST=''; OUTPUT=''
THREADS=8; MINREADS=50; MINTPM=1.0; SALMON_INDEX_OVERRIDE=''

OPTS=$(getopt -o 1:2:o:p: \
  -l input1:,input2:,host:,output:,threads:,minreads:,mintpm:,salmon-index: \
  -- "$@")
eval set -- "$OPTS"

while [ $# -gt 0 ]; do
  case "$1" in
    -1|--input1)        INPUT1=$2; shift;;
    -2|--input2)        INPUT2=$2; shift;;
    --host)             HOST=$2; shift;;
    -o|--output)        OUTPUT=$2; shift;;
    -p|--threads)       THREADS=$2; shift;;
    --minreads)         MINREADS=$2; shift;;
    --mintpm)           MINTPM=$2; shift;;
    --salmon-index)     SALMON_INDEX_OVERRIDE=$2; shift;;
    (--)                shift; break;;
  esac
  shift
done

#  Validation
if [ -z "$INPUT1" ] || [ -z "$HOST" ] || [ -z "$OUTPUT" ]; then
  echo "ERROR: Missing required parameters (-1, --host, -o)"; exit 1
fi
[ -f "$INPUT1" ] || { echo "ERROR: Input file not found: $INPUT1"; exit 1; }

if [ -z "$INPUT2" ]; then
  MODE="single-end"
  echo "Running in SINGLE-END mode"
else
  [ -f "$INPUT2" ] || { echo "ERROR: Input file 2 not found: $INPUT2"; exit 1; }
  MODE="paired-end"
  echo "Running in PAIRED-END mode"
fi

#  Paths 
BASEDIR="$(cd "$(dirname "$0")" && pwd)"
DBDIR="$BASEDIR/database/bowtie2"

case $HOST in
  cattle)
    HOST_INDEX="bos_taurus_index"
    MICROBIOME_INDEX="rumen_microbiome_index"
    FUNGAL_INDEX=""
    VIRAL_SCREEN="all_viral_index"
    ;;
  pig)
    HOST_INDEX="pig_genome_index"
    MICROBIOME_INDEX="pig_microbiome_index"
    FUNGAL_INDEX=""
    VIRAL_SCREEN="all_viral_index"
    ;;
  grapevine)
    HOST_INDEX="grapevine_genome_index"
    MICROBIOME_INDEX="grapevine_microbiome_index"
    FUNGAL_INDEX=""
    VIRAL_SCREEN="all_viral_index"
    ;;
  human)
    HOST_INDEX="human_genome_index"
    MICROBIOME_INDEX="hrgmv2_hybrid_bt2"
    FUNGAL_INDEX="cgf_top100_bt2"
    VIRAL_SCREEN="all_viral_index"
    ;;
  *)
    echo "ERROR: Invalid host. Choose: cattle, pig, grapevine, human"; exit 1;;
esac

# Salmon index path — human uses curated index, other hosts use full viral DB
if [ -n "$SALMON_INDEX_OVERRIDE" ]; then
SALMON_IDX="$SALMON_INDEX_OVERRIDE"
elif [ "$HOST" = "human" ]; then
SALMON_IDX="$BASEDIR/database/salmon_index"
else
SALMON_IDX="$BASEDIR/database/salmon_index_all"
fi
[ -d "$SALMON_IDX" ] || { echo "ERROR: Salmon index not found: $SALMON_IDX"; exit 1; }

mkdir -p "$OUTPUT"
echo "Host: $HOST | Mode: $MODE | Threads: $THREADS"
echo ""

# STEP 1 - Broad viral pre-screen
echo "Step 1: Broad viral pre-screen..."
if [ "$MODE" = "paired-end" ]; then
  bowtie2 -x "$DBDIR/$VIRAL_SCREEN" \
    -1 "$INPUT1" -2 "$INPUT2" \
    --sensitive-local --no-unal \
    -S "$OUTPUT/step1_viral.sam" -p "$THREADS"
  samtools view -b "$OUTPUT/step1_viral.sam" | \
    samtools fastq -N \
      -1 "$OUTPUT/step1_R1.fastq" \
      -2 "$OUTPUT/step1_R2.fastq" \
      -0 /dev/null -s /dev/null -n -
  STEP1_R1="$OUTPUT/step1_R1.fastq"
  STEP1_R2="$OUTPUT/step1_R2.fastq"
else
  bowtie2 -x "$DBDIR/$VIRAL_SCREEN" \
    -U "$INPUT1" --sensitive-local --no-unal \
    -S "$OUTPUT/step1_viral.sam" -p "$THREADS"
  grep -v ^@ "$OUTPUT/step1_viral.sam" | \
    awk '{print "@"$1"\n"$10"\n+\n"$11}' > "$OUTPUT/step1_viral.fastq"
  STEP1_R1="$OUTPUT/step1_viral.fastq"
  STEP1_R2=""
fi
rm -f "$OUTPUT/step1_viral.sam"

# STEP 2 - Host genome filter
echo "Step 2: Filtering host genome..."
if [ "$MODE" = "paired-end" ]; then
  bowtie2 -x "$DBDIR/$HOST_INDEX" \
    -1 "$STEP1_R1" -2 "$STEP1_R2" \
    --un-conc "$OUTPUT/step2_R%.fastq" \
    -p "$THREADS" -S /dev/null
  STEP2_R1="$OUTPUT/step2_R1.fastq"
  STEP2_R2="$OUTPUT/step2_R2.fastq"
else
  bowtie2 -x "$DBDIR/$HOST_INDEX" \
    -U "$STEP1_R1" --un "$OUTPUT/step2.fastq" \
    -p "$THREADS" -S /dev/null
  STEP2_R1="$OUTPUT/step2.fastq"
  STEP2_R2=""
fi

# STEP 3 - Microbiome filter
echo "Step 3: Filtering microbiome..."
if [ "$MODE" = "paired-end" ]; then
  bowtie2 -x "$DBDIR/$MICROBIOME_INDEX" \
    -1 "$STEP2_R1" -2 "$STEP2_R2" \
    --un-conc "$OUTPUT/step3_R%.fastq" \
    -p "$THREADS" -S /dev/null
  STEP3_R1="$OUTPUT/step3_R1.fastq"
  STEP3_R2="$OUTPUT/step3_R2.fastq"
else
  bowtie2 -x "$DBDIR/$MICROBIOME_INDEX" \
    -U "$STEP2_R1" --un "$OUTPUT/step3.fastq" \
    -p "$THREADS" -S /dev/null
  STEP3_R1="$OUTPUT/step3.fastq"
  STEP3_R2=""
fi

# STEP 4 - Fungal filter (human only)
if [ -n "$FUNGAL_INDEX" ]; then
  echo "Step 4: Filtering fungi..."
  if [ "$MODE" = "paired-end" ]; then
    bowtie2 -x "$DBDIR/$FUNGAL_INDEX" \
      -1 "$STEP3_R1" -2 "$STEP3_R2" \
      --un-conc "$OUTPUT/step4_R%.fastq" \
      -p "$THREADS" -S /dev/null
    FINAL_R1="$OUTPUT/step4_R1.fastq"
    FINAL_R2="$OUTPUT/step4_R2.fastq"
  else
    bowtie2 -x "$DBDIR/$FUNGAL_INDEX" \
      -U "$STEP3_R1" --un "$OUTPUT/step4.fastq" \
      -p "$THREADS" -S /dev/null
    FINAL_R1="$OUTPUT/step4.fastq"
    FINAL_R2=""
  fi
else
  echo "Step 4: Skipping fungal filter (non-human host)"
  FINAL_R1="$STEP3_R1"
  FINAL_R2="$STEP3_R2"
fi

# STEP 5 - Salmon quantification
echo "Step 5: Salmon EM quantification..."
mkdir -p "$OUTPUT/salmon_quant"

if [ "$MODE" = "paired-end" ]; then
  salmon quant \
    -i "$SALMON_IDX" \
    -l A \
    -1 "$FINAL_R1" \
    -2 "$FINAL_R2" \
    --validateMappings \
    --rangeFactorizationBins 4 \
    --minAssignedFrags 1 \
    -p "$THREADS" \
    -o "$OUTPUT/salmon_quant" \
    --quiet
else
  salmon quant \
    -i "$SALMON_IDX" \
    -l A \
    -r "$FINAL_R1" \
    --validateMappings \
    --rangeFactorizationBins 4 \
    --minAssignedFrags 1 \
    -p "$THREADS" \
    -o "$OUTPUT/salmon_quant" \
    --quiet
fi

# STEP 6 - Filter quant.sf and generate report
echo "Step 6: Filtering detections..."
QUANT="$OUTPUT/salmon_quant/quant.sf"

# Filtered detections
awk -v minr="$MINREADS" -v mint="$MINTPM" \
  'NR>1 && ($5+0)>=minr && ($4+0)>=mint {print $1"\t"$5"\t"$4}' \
  "$QUANT" | sort -k2 -rn > "$OUTPUT/detected_viruses.txt"

# Full unfiltered output
awk 'NR>1 && ($5+0)>0 {print $1"\t"$5"\t"$4}' \
  "$QUANT" | sort -k2 -rn > "$OUTPUT/detected_viruses_all.txt"

# Human blacklist (EVEs, integrated elements)
if [ "$HOST" = "human" ]; then
  BLACKLIST="NC_022518.1|NC_021858.1|NC_006658.1|NC_032111.1|AC_000005.1|AC_000006.1|AC_000007.1|AC_000008.1|AC_000017.1|AC_000018.1|AC_000019.1|AC_000010.1|AC_000011.1"
  grep -vE "$BLACKLIST" "$OUTPUT/detected_viruses.txt" \
    > "$OUTPUT/detected_viruses_tmp.txt" || true
  mv "$OUTPUT/detected_viruses_tmp.txt" "$OUTPUT/detected_viruses.txt"
fi

N_VIRUSES=$(wc -l < "$OUTPUT/detected_viruses.txt")
TOTAL_READS=$(( $(wc -l < "$INPUT1") / 4 ))
MAPPED=$(awk 'NR>1 {sum+=$5} END {print int(sum)}' "$QUANT")

cat > "$OUTPUT/REPORT.txt" << EOF
viromescan v2.2 Analysis Report
================================
Mode:              $MODE
Host:              $HOST
Input R1:          $INPUT1
Total input reads: $TOTAL_READS
Salmon index:      $SALMON_IDX

Parameters:
  Threads:         $THREADS
  Min NumReads:    $MINREADS
  Min TPM:         $MINTPM

Results:
  Total viral reads assigned: $MAPPED
  Viruses detected:           $N_VIRUSES

Top detected viruses (Name / NumReads / TPM):
$(head -20 "$OUTPUT/detected_viruses.txt" | \
  awk '{printf "  %-30s reads=%-10s TPM=%s\n", $1, $2, $3}')

Output files:
  $OUTPUT/detected_viruses.txt      — filtered detections
  $OUTPUT/detected_viruses_all.txt  — all with >0 reads
  $OUTPUT/salmon_quant/quant.sf     — full Salmon output
  $OUTPUT/REPORT.txt                — this report

Analysis completed: $(date)
EOF

# Cleanup intermediates
rm -f "$OUTPUT"/step*.fastq "$OUTPUT"/*.sam

echo ""
echo "viromescan v2.2 complete. Results in: $OUTPUT/"
echo "Viruses detected: $N_VIRUSES"
echo "Full report: $OUTPUT/REPORT.txt"
