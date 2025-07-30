#!/usr/bin/env bash
# AraFusion.sh
# “A Complete Hybrid Assembly Workflow”
# Usage: bash AraFusion.sh
# Author & Contact : ykaya@mpipz.mpg.de | yyasinkkaya@gmail.com

set -euo pipefail
IFS=$'\n\t'

#### 0) TOOLS & PARAMETERS ##########################################

# Hybrid assembly
VERKKO="/home/ykaya/.conda/envs/verkko/bin/verkko"
HIFIASM="/netscratch/irg/grp_hancock/Raw_HiFi_Data/mpgc_Arabidopsis_EA/version2/hifiasm/hifiasm"
TELO_M="TTTAGGG"

# Downstream tools
NUCMER=/opt/share/software/bin/nucmer
DELTA_FILTER=/opt/share/software/bin/delta-filter
QUICKMERGE=/opt/share/software/bin/quickmerge
PBCSTAT=pbcstat
CALCUTS=calcuts
HIST_PLOT=hist_plot.py
SPLIT_FA=split_fa
MINIMAP2=minimap2
PURGE_DUPS=purge_dups
GET_SEQS=get_seqs
RAGTAG=ragtag.py
SAMTOOLS=samtools
MAKE_AGP=make_agp.py
CLOSE_GAPS=/xxx/MaSuRCA-4.1.0/bin/close_scaffold_gaps.sh
PILON=pilon

# Threads & parameters
THREADS_HIFI=50     # Verkko2
THREADS_HIA=32      # hifiasm
THREADS=20          # downstream steps
FILTER_LEN=10000
MIN_ID=0.9
HCO=5.0
C_VAL=1.5
MIN_ALIGN=5000

# Directories
HIFI_READS="/path/to/HiFi_reads"
ONT_READS="/path/to/ONT_reads"
COLCC_REF="/path/to/Col-CC_reference.fasta"
FIRST_OUT="/path/to/genome_assembly/first_assembly"

ASSEMBLY_DIR="/path/to/genome_assembly"
NGS_DIR="/path/to/Illumina_allData"

# Samples with both HiFi & ONT for initial hybrid assembly
SAMPLES_WITH_ONT=( ETx,  ETy, ETz, ...)


#### 1) HYBRID ASSEMBLY WITH VERKKO2 #################################
mkdir -p "$FIRST_OUT/verkko2"
for sample in "${SAMPLES_WITH_ONT[@]}"; do
  echo "=== Verkko2 hybrid assembly: $sample ==="
  OUTDIR="$FIRST_OUT/verkko2/$sample"
  mkdir -p "$OUTDIR"
  $VERKKO -d "$OUTDIR" \
    --hifi "$HIFI_READS/${sample}.hf.fastq.gz" \
    --nano "$ONT_READS/${sample}.ONT.fastq.gz" \
    --ref "$COLCC_REF" \
    --telomere-motif "$TELO_M" \
    --grid --local-memory 400 --local-cpus $THREADS_HIFI
  echo "--- Completed Verkko2 for $sample"
done

#### 2) HYBRID ASSEMBLY WITH HIFIASM ##################################
mkdir -p "$FIRST_OUT/hifiasm"
for sample in "${SAMPLES_WITH_ONT[@]}"; do
  echo "=== hifiasm hybrid assembly: $sample ==="
  OUT_PREFIX="$FIRST_OUT/hifiasm/${sample}.asm"
  HIFI_FASTQ=$(ls "$HIFI_READS/${sample}".fastq.gz "$HIFI_READS/${sample}".fq.gz 2>/dev/null | head -n1)
  ONT_FASTQ=$(ls "$ONT_READS/${sample}".fastq "$ONT_READS/${sample}".fastq.gz 2>/dev/null | head -n1)
  [[ -z "$HIFI_FASTQ" ]] && echo "WARN: missing HiFi for $sample" && continue
  [[ -z "$ONT_FASTQ" ]] && echo "WARN: missing ONT for $sample" && continue

  $HIFIASM \
    -t $THREADS_HIA -o "$OUT_PREFIX" \
    -l0 --hg-size 140m \
    -ul "$ONT_FASTQ" --ul-cut 100000 \
    --dual-scaf --primary --telo-m "$TELO_M" \
    "$HIFI_FASTQ"
  echo "--- Completed hifiasm for $sample"
done


#Samples for merging & downstream (define after assemblies to allow N50 calculation)
SAMPLES=( ETx1 ETx2 ETx3 )
declare -A N50_CUTOFF=( [ETx1]=12000000 [ETx2]=9605000 [ETx3]=16039000 )


#### 3) NUCMER ALIGNMENT & QUICKMERGE ###############################
for sample in "${SAMPLES[@]}"; do
  echo "=== QuickMerge: $sample ==="
  BASE="$ASSEMBLY_DIR/$sample"
  ASM="$FIRST_OUT/hifiasm/${sample}.asm.p_ctg.fa"
  VERK2="$FIRST_OUT/verkko2/$sample/assembly.fasta"
  QM_DIR="$BASE/05_quickmerge"
  PREFIX="${sample}.asm_verkko"
  mkdir -p "$QM_DIR"

  # 3.1) alignment
  $NUCMER -t $THREADS -p "$QM_DIR/$PREFIX.aln" "$ASM" "$VERK2"
  $DELTA_FILTER -r -q -l $FILTER_LEN -i $MIN_ID \
    "$QM_DIR/$PREFIX.aln.delta" > "$QM_DIR/$PREFIX.aln.flt.delta"

  # 3.2) merge
  $QUICKMERGE \
    -d "$QM_DIR/$PREFIX.aln.flt.delta" \
    -q "$VERK2" \
    -r "$ASM" \
    -hco $HCO -c $C_VAL \
    -l ${N50_CUTOFF[$sample]} -ml $MIN_ALIGN \
    -p "$PREFIX.merge" -v > "$QM_DIR/quickmerge.log" 2>&1

  MERGED="$QM_DIR/$PREFIX.merge.fasta"

  #### 4) PURGE HAPLOTIGS ##########################################
  echo "--- Purging haplotigs: $sample"
  PDIR="$QM_DIR/purge_dups"; mkdir -p "$PDIR"
  $PBCSTAT -O "$PDIR" "$MERGED.nano.paf"
  $CALCUTS "$PDIR/PB.stat" > "$PDIR/cutoffs.auto"
  $HIST_PLOT -X 250 -c "$PDIR/cutoffs.auto" \
    "$PDIR/PB.stat" "$PDIR/PB.stat.auto.png"
  $SPLIT_FA "$MERGED" > "$PDIR/${sample}.split.fa"
  $MINIMAP2 -x asm5 -DP -t $THREADS \
    "$PDIR/${sample}.split.fa" "$PDIR/${sample}.split.fa" \
    -o "$PDIR/${sample}.self.paf"
  $PURGE_DUPS -2 -T "$PDIR/cutoffs.auto" -c "$PDIR/PB.base.cov" \
    "$PDIR/${sample}.self.paf" > "$PDIR/dups.bed"
  $GET_SEQS -e -p "${sample}.purged" "$PDIR/dups.bed" "$MERGED"

  PURGED_FA="$PDIR/${sample}.purged.fa"

  #### 5) SCAFFOLDING ##############################################
  echo "--- RagTag scaffold: $sample"
  SDIR="$BASE/06_scaffolding"; mkdir -p "$SDIR" && cd "$SDIR"
  ln -sf "$PURGED_FA" .; ln -sf "$COLCC_REF" reference.fa
  $RAGTAG scaffold reference.fa $(basename "$PURGED_FA") -t $THREADS > ragtag.log 2>&1
  cd "$BASE"

  #### 6) GAP-CLOSING #############################################
  echo "--- Gap-closing (HiFi): $sample"
  GDIR="$BASE/07_gapclosing"; mkdir -p "$GDIR" && cd "$GDIR"
  ln -sf "$FIRST_OUT/verkko2/$sample/assembly.fasta" asm.fa
  ln -sf "$FIRST_OUT/verkko2/$sample/hifi.corrected.fasta.gz" . && gunzip hifi.corrected.fasta.gz
  cd gap_closing_reads
  $CLOSE_GAPS -r ../asm.fa -q ../hifi.corrected.fasta -t $THREADS

  #### 7) POLISH with Pilon #######################################
  echo "--- Pilon polishing: $sample"
  PDIR2="$BASE/08_pilon"; mkdir -p "$PDIR2" && cd "$PDIR2"
  if [[ -f "$NGS_DIR/${sample}_R1.fastq.gz" ]]; then
    R1="$NGS_DIR/${sample}_R1.fastq.gz"; R2="$NGS_DIR/${sample}_R2.fastq.gz"
  elif [[ -f "$NGS_DIR/${sample}_R1.fastq" ]]; then
    R1="$NGS_DIR/${sample}_R1.fastq";   R2="$NGS_DIR/${sample}_R2.fastq"
  else
    echo "--- Missing NGS for $sample; skipping Pilon"; cd "$BASE"; continue
  fi
  inputAsm="$GDIR/gap_closing_reads/asm.closed.fa"
  for r in {1..3}; do
    SAM="$PDIR2/${sample}.r${r}.sam"
    BAM="$PDIR2/${sample}.r${r}.sorted.bam"
    $MINIMAP2 -ax sr -t $THREADS "$inputAsm" "$R1" "$R2" > "$SAM"
    $SAMTOOLS view -Sb "$SAM" | $SAMTOOLS sort -o "$BAM"
    $SAMTOOLS index "$BAM"
    $PILON --genome "$inputAsm" --frags "$BAM" --outdir "$PDIR2" --output "${sample}.pilon${r}"
    inputAsm="$PDIR2/${sample}.pilon${r}.fasta"
  done
  echo "--- Completed $sample: $inputAsm"
  cd "$BASE"
done

echo "ALL SAMPLES COMPLETE 🎉"