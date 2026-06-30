#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script dir is: $SCRIPT_DIR"

usage() {
  cat <<EOF
Usage:
  $0 -g GENOME_DIR -s SAMPLE_LIST -a ALIGNER -c CANCER_TYPE [-i INPUT_DIR] [-o OUTPUT_DIR] [-t THREADS]

Required:
  -g GENOME_DIR   Path to genome folder containing:
                  GRCh38.primary_assembly.genome.fa
                  gencode.v49.annotation.gtf
                  dbsnp_146.hg38.vcf.gz

                  For STAR mode:
                  STAR/   (STAR genome index directory)

                  For HISAT2 mode:
                  GRCh38.primary_assembly.genome.*.ht2

  -s SAMPLE_LIST     File with one sample ID per line

  -a ALIGNER      Aligner to use: star or hisat2

  -c CANCER_TYPE  breast or ovary

Optional:
  -i INPUT_DIR    Directory containing input FASTQ files
                  Default: current directory

  -o OUTPUT_DIR   Directory for results
                  Default: ./RNA_variant_pipeline_out

  -t THREADS      Number of threads
                  Default: \$SLURM_CPUS_PER_TASK if set, otherwise 8

Example:
$0 -g /path/to/genome -s sample_ids.txt -a star   -c breast -i /path/to/fastq -o results -t 16
$0 -g /path/to/genome -s sample_ids.txt -a star   -c ovary  -i /path/to/fastq -o results -t 16

$0 -g /path/to/genome -s sample_ids.txt -a hisat2 -c breast -i /path/to/fastq -o results -t 16
$0 -g /path/to/genome -s sample_ids.txt -a hisat2 -c ovary  -i /path/to/fastq -o results -t 16
EOF
  exit 1
}

# defaults

INPUT_DIR="."
OUTPUT_DIR="./RNA_variant_pipeline_out"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
ALIGNER=""
CANCER_TYPE=""

while getopts ":g:s:a:c:i:o:t:h" opt; do
  case $opt in
    g) GENOME_DIR="$OPTARG" ;;
    s) SAMPLE_LIST="$OPTARG" ;;
    a) ALIGNER="$OPTARG" ;;
    i) INPUT_DIR="$OPTARG" ;;
    c) CANCER_TYPE="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage ;;
    \?) echo "Error: Invalid option -$OPTARG" >&2; usage ;;
    :) echo "Error: Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done


# required args

[[ -z "${GENOME_DIR:-}" ]] && { echo "Error: -g GENOME_DIR is required"; usage; }
[[ -z "${SAMPLE_LIST:-}"   ]] && { echo "Error: -s SAMPLE_LIST is required"; usage; }
[[ -z "${ALIGNER:-}"    ]] && { echo "Error: -a ALIGNER is required"; usage; }

if [[ "$ALIGNER" != "star" && "$ALIGNER" != "hisat2" ]]; then
  echo "Error: -a must be either 'star' or 'hisat2'"
  exit 1
fi
[[ -z "${CANCER_TYPE:-}" ]] && { echo "Error: -c CANCER_TYPE is required"; usage; }

if [[ "$CANCER_TYPE" != "breast" && "$CANCER_TYPE" != "ovary" ]]; then
  echo "Error: -c must be either 'breast' or 'ovary'"
  exit 1
fi


# genome files

REF="${GENOME_DIR}/GRCh38.primary_assembly.genome.fa"
GTF="${GENOME_DIR}/gencode.v49.annotation.gtf"
DBSNP="${GENOME_DIR}/dbsnp_146.hg38.vcf.gz"

STAR_INDEX="${GENOME_DIR}/STAR"
HISAT_INDEX="${GENOME_DIR}/GRCh38.primary_assembly.genome"


# checks

[[ -d "$GENOME_DIR" ]] || { echo "Error: Genome directory not found: $GENOME_DIR"; exit 1; }
[[ -f "$SAMPLE_LIST" ]]   || { echo "Error: Sample list file not found: $SAMPLE_LIST"; exit 1; }
[[ -f "$REF" ]]        || { echo "Error: Reference fasta not found: $REF"; exit 1; }
[[ -f "$GTF" ]]        || { echo "Error: GTF not found: $GTF"; exit 1; }
[[ -f "$DBSNP" ]]      || { echo "Error: dbSNP VCF not found: $DBSNP"; exit 1; }

if [[ "$ALIGNER" == "star" ]]; then
  [[ -d "$STAR_INDEX" ]] || { echo "Error: STAR index directory not found: $STAR_INDEX"; exit 1; }
fi

if [[ "$ALIGNER" == "hisat2" ]]; then
  ls "${HISAT_INDEX}".*.ht2 >/dev/null 2>&1 || {
    echo "Error: HISAT2 index files not found with basename: $HISAT_INDEX"
    exit 1
  }
fi


# tool checks

COMMON_TOOLS=(fastp samtools stringtie bgzip tabix bcftools python vcf-sort)

for cmd in "${COMMON_TOOLS[@]}"; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "Error: Required command not found in PATH: $cmd"
    exit 1
  }
done

if [[ "$ALIGNER" == "star" ]]; then
  command -v STAR >/dev/null 2>&1 || {
    echo "Error: Required command not found in PATH: STAR"
    exit 1
  }
fi

if [[ "$ALIGNER" == "hisat2" ]]; then
  command -v hisat2 >/dev/null 2>&1 || {
    echo "Error: Required command not found in PATH: hisat2"
    exit 1
  }
fi


# helper functions

file_ok() {
  [[ -s "$1" ]]
}

log_step() {
  echo "[$(date '+%F %T')] $1"
}

bam_ok() {
  local bam="$1"
  [[ -s "$bam" ]] || return 1
  samtools quickcheck -q "$bam" || return 1
  return 0
}


find_fastq() {
  local base="$1"
  for ext in fastq fastq.gz fq fq.gz; do
    if [[ -f "${INPUT_DIR}/${base}.${ext}" ]]; then
      echo "${INPUT_DIR}/${base}.${ext}"
      return 0
    fi
  done
  return 1
}

fastq_ok() {
  local f="$1"

  [[ -s "$f" ]] || return 1

  # gzip must be valid
  gzip -t "$f" 2>/dev/null || return 1

  # first record header must start with @
  zcat "$f" | head -n 1 | grep -q '^@' || return 1

  # first 4000 lines must look like complete FASTQ records
  zcat "$f" | head -n 4000 | awk '
    NR%4==1 { if ($0 !~ /^@/) exit 1 }
    NR%4==3 { if ($0 !~ /^\+/) exit 1 }
    END { if (NR<4 || NR%4!=0) exit 1 }
  ' || return 1

  return 0
}


# output dirs

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/ST_gtf"
mkdir -p "$OUTPUT_DIR/logs"

echo "=========================================="
echo "Genome dir : $GENOME_DIR"
echo "Sample list   : $SAMPLE_LIST"
echo "Input dir  : $INPUT_DIR"
echo "Output dir : $OUTPUT_DIR"
echo "Threads    : $THREADS"
echo "Aligner    : $ALIGNER"
echo "Cancer type    : $CANCER_TYPE"
if [[ "$ALIGNER" == "star" ]]; then
  echo "STAR index : $STAR_INDEX"
else
  echo "HISAT2 idx : $HISAT_INDEX"
fi
echo "=========================================="


# main loop

while read -r srr || [[ -n "$srr" ]]; do
  [[ -z "$srr" ]] && continue

  log_step "Processing $srr"

  R1=$(find_fastq "${srr}_1" || true)
  R2=$(find_fastq "${srr}_2" || true)

  if [[ -z "$R1" || -z "$R2" ]]; then
    log_step "Warning: FASTQ files missing for $srr, skipping"
    continue
  fi

  SAMPLE_DIR="${OUTPUT_DIR}/${srr}"
  mkdir -p "$SAMPLE_DIR"
  mkdir -p "${SAMPLE_DIR}/ballgown_table/${srr}"
  mkdir -p "${OUTPUT_DIR}/ST_gtf/${srr}"

  FINAL_OUT="${SAMPLE_DIR}/${srr}_annotated.vcf"

if file_ok "$FINAL_OUT"; then
  log_step "$srr: final annotated VCF exists, skipping upstream steps"

  # still run SNP extraction
  SNP_OUT="${SAMPLE_DIR}/${srr}_${CANCER_TYPE}_selected_SNP.txt"

  if [[ "$CANCER_TYPE" == "breast" ]]; then
    SNP_PATTERN="rs2366152|rs7091441"
  else
    SNP_PATTERN="rs932501|rs2072588|rs7318592|rs9506960|rs9510420|rs12583808|rs60135126"
  fi

  if file_ok "$SNP_OUT"; then
    log_step "$srr: tissue-specific SNP list already extracted, skipping"
  else
    log_step "$srr: extracting $CANCER_TYPE tissue-specific SNPs"
    egrep -w "$SNP_PATTERN" "$FINAL_OUT" | cut -f3 > "$SNP_OUT" || true
  fi

  continue
fi


# 1. fastp
if bam_ok "${SAMPLE_DIR}/${srr}_sorted.bam" && file_ok "${SAMPLE_DIR}/${srr}_sorted.bam.bai"; then
  log_step "$srr: sorted BAM + index already exist, skipping fastp"
elif fastq_ok "${SAMPLE_DIR}/${srr}_1_trimmed.fastq.gz" && fastq_ok "${SAMPLE_DIR}/${srr}_2_trimmed.fastq.gz"; then
  log_step "$srr: fastp already done, skipping"
else
  log_step "$srr: running fastp"

  rm -f "${SAMPLE_DIR}/${srr}_1_trimmed.fastq.gz" \
        "${SAMPLE_DIR}/${srr}_2_trimmed.fastq.gz" \
        "${SAMPLE_DIR}/${srr}.fastp.json" \
        "${SAMPLE_DIR}/${srr}.fastp.html"

  fastp \
    -i "$R1" -I "$R2" \
    -o "${SAMPLE_DIR}/${srr}_1_trimmed.fastq.gz" \
    -O "${SAMPLE_DIR}/${srr}_2_trimmed.fastq.gz" \
    --cut_front --cut_tail \
    --cut_mean_quality 30 \
    --length_required 25 \
    --thread "$THREADS" \
    --json "${SAMPLE_DIR}/${srr}.fastp.json" \
    --detect_adapter_for_pe \
    --html "${SAMPLE_DIR}/${srr}.fastp.html"
fi



  # 2 + 3. Alignment + sorted BAM + index
  if file_ok "${SAMPLE_DIR}/${srr}_sorted.bam" && file_ok "${SAMPLE_DIR}/${srr}_sorted.bam.bai"; then
    log_step "$srr: sorted BAM already done, skipping alignment"
  else
    if [[ "$ALIGNER" == "star" ]]; then
      log_step "$srr: running STAR"
      STAR \
        --genomeDir "$STAR_INDEX" \
        --runThreadN "$THREADS" \
        --readFilesCommand gunzip -c \
        --readFilesIn "${SAMPLE_DIR}/${srr}_1_trimmed.fastq.gz" "${SAMPLE_DIR}/${srr}_2_trimmed.fastq.gz" \
        --outSAMattributes NH HI AS nM XS \
        --outFileNamePrefix "${SAMPLE_DIR}/${srr}_" \
        --outSAMtype BAM SortedByCoordinate \
        > "${OUTPUT_DIR}/logs/${srr}.STAR.log" 2>&1

      log_step "$srr: running samtools calmd"
      samtools calmd -b \
        "${SAMPLE_DIR}/${srr}_Aligned.sortedByCoord.out.bam" \
        "$REF" \
        > "${SAMPLE_DIR}/${srr}_MD.bam"

      log_step "$srr: sorting final BAM"
      samtools sort -@ "$THREADS" \
        "${SAMPLE_DIR}/${srr}_MD.bam" \
        -o "${SAMPLE_DIR}/${srr}_sorted.bam"

      samtools index "${SAMPLE_DIR}/${srr}_sorted.bam"

      rm -f "${SAMPLE_DIR}/${srr}_Aligned.sortedByCoord.out.bam"
      rm -f "${SAMPLE_DIR}/${srr}_MD.bam"

    elif [[ "$ALIGNER" == "hisat2" ]]; then
      log_step "$srr: running HISAT2"
      hisat2 -p "$THREADS" -q \
        -x "$HISAT_INDEX" \
        -1 "${SAMPLE_DIR}/${srr}_1_trimmed.fastq.gz" \
        -2 "${SAMPLE_DIR}/${srr}_2_trimmed.fastq.gz" \
        -S "${SAMPLE_DIR}/${srr}.sam" \
        2> "${OUTPUT_DIR}/logs/${srr}.hisat2.log"

      log_step "$srr: converting/sorting BAM"
      samtools view -bS -@ "$THREADS" "${SAMPLE_DIR}/${srr}.sam" | \
        samtools sort -@ "$THREADS" -o "${SAMPLE_DIR}/${srr}_sorted.bam"

      samtools index "${SAMPLE_DIR}/${srr}_sorted.bam"
      rm -f "${SAMPLE_DIR}/${srr}.sam"
    fi
  fi

  # 4. StringTie
  if file_ok "${OUTPUT_DIR}/ST_gtf/${srr}/${srr}.gtf"; then
    log_step "$srr: StringTie already done, skipping"
  else
    log_step "$srr: running StringTie"
    stringtie \
      -b "${SAMPLE_DIR}/ballgown_table/${srr}/" \
      -p "$THREADS" \
      -G "$GTF" \
      -o "${OUTPUT_DIR}/ST_gtf/${srr}/${srr}.gtf" \
      "${SAMPLE_DIR}/${srr}_sorted.bam"
  fi

  # 5. Extract FPKM
  if file_ok "${SAMPLE_DIR}/${srr}_FPKM.txt"; then
    log_step "$srr: FPKM already extracted, skipping"
  else
    if [[ -f "${SAMPLE_DIR}/ballgown_table/${srr}/t_data.ctab" ]]; then
      log_step "$srr: extracting FPKM"
      cut -f6,12 "${SAMPLE_DIR}/ballgown_table/${srr}/t_data.ctab" > "${SAMPLE_DIR}/${srr}_FPKM.txt"
    else
      log_step "$srr: t_data.ctab not found after StringTie, skipping FPKM extraction"
    fi
  fi

  # 6. Opossum
if bam_ok "${SAMPLE_DIR}/${srr}_opossum.bam" && file_ok "${SAMPLE_DIR}/${srr}_opossum.bam.bai"; then
    log_step "$srr: Opossum already done, skipping"
  else
    if [[ "$ALIGNER" == "star" ]]; then
      log_step "$srr: running Opossum for STAR (SoftClipsExist=True)"
      echo ">>> INFO: STAR detected → running Opossum with --SoftClipsExist True"

      python "${SCRIPT_DIR}/Opossum/Opossum.py" \
        --BamFile "${SAMPLE_DIR}/${srr}_sorted.bam" \
        --OutFile "${SAMPLE_DIR}/${srr}_opossum.bam" \
        --SoftClipsExist True

    elif [[ "$ALIGNER" == "hisat2" ]]; then
      log_step "$srr: running Opossum for HISAT2 (default settings)"
      echo ">>> INFO: HISAT2 detected → running Opossum without soft clip flag"

      python "${SCRIPT_DIR}/Opossum/Opossum.py" \
        --BamFile "${SAMPLE_DIR}/${srr}_sorted.bam" \
        --OutFile "${SAMPLE_DIR}/${srr}_opossum.bam"
    fi
  fi


  # 7. Platypus
  if file_ok "${SAMPLE_DIR}/${srr}_platypus_variants.vcf"; then
    log_step "$srr: Platypus already done, skipping"
  else
    log_step "$srr: running Platypus"
    python "${SCRIPT_DIR}/Platypus/bin/Platypus.py" callVariants \
      --bamFiles "${SAMPLE_DIR}/${srr}_opossum.bam" \
      --refFile "$REF" \
      --filterDuplicates 0 \
      --minReads 5 \
      --minMapQual 20 \
      --minFlank 0 \
      --maxReadLength 500 \
      --minGoodQualBases 20 \
      --minBaseQual 30 \
      --verbosity 3 \
      -o "${SAMPLE_DIR}/${srr}_platypus_variants.vcf"
  fi

  # 8 + 9. sort/compress/index VCF
  if file_ok "${SAMPLE_DIR}/${srr}_sorted.vcf.gz" && file_ok "${SAMPLE_DIR}/${srr}_sorted.vcf.gz.tbi"; then
    log_step "$srr: sorted/indexed VCF already done, skipping"
  else
    log_step "$srr: sorting/compressing/indexing VCF"
    vcf-sort "${SAMPLE_DIR}/${srr}_platypus_variants.vcf" > "${SAMPLE_DIR}/${srr}_sorted.vcf"
    bgzip -c "${SAMPLE_DIR}/${srr}_sorted.vcf" > "${SAMPLE_DIR}/${srr}_sorted.vcf.gz"
    tabix -p vcf "${SAMPLE_DIR}/${srr}_sorted.vcf.gz"
  fi

  # 10. annotate with dbSNP
  if file_ok "${SAMPLE_DIR}/${srr}_annotated.vcf"; then
    log_step "$srr: annotation already done, skipping"
  else
    log_step "$srr: annotating with dbSNP"
    bcftools annotate \
      -c CHROM,FROM,TO,ID,INFO \
      -a "$DBSNP" \
      -o "${SAMPLE_DIR}/${srr}_annotated.vcf" \
      "${SAMPLE_DIR}/${srr}_sorted.vcf.gz"
  fi

  # 11. select tissue-specific SNP IDs
  SNP_OUT="${SAMPLE_DIR}/${srr}_${CANCER_TYPE}_selected_SNP.txt"

  if [[ "$CANCER_TYPE" == "breast" ]]; then
    SNP_PATTERN="rs2366152|rs7091441"
  else
    SNP_PATTERN="rs932501|rs2072588|rs7318592|rs9506960|rs9510420|rs12583808|rs60135126"
  fi

  if file_ok "$SNP_OUT"; then
    log_step "$srr: tissue-specific SNP list already extracted, skipping"
  else
    log_step "$srr: extracting $CANCER_TYPE tissue-specific SNPs"
    egrep -w "$SNP_PATTERN" "${SAMPLE_DIR}/${srr}_annotated.vcf" | cut -f3 > "$SNP_OUT" || true
  fi


  log_step "$srr DONE ✅"

done < "$SAMPLE_LIST"
