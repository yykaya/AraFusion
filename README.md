# AraFusion: A Hybrid Genome Assembly Pipeline for Arabidopsis thaliana

## Introduction
Advances in long-read sequencing have made telomere-to-telomere (T2T) assemblies increasingly routine, as technology costs decline and bioinformatics tools become ever more versatile. Integrating different sequencing data types—Oxford Nanopore (ONT) ultra-long reads, PacBio HiFi, and paired-end short-read polishing—enables both high contiguity and base-level accuracy. **AraFusion** unites these approaches into a streamlined workflow tailored for _Arabidopsis thaliana_, delivering chromosome-scale, high-fidelity genomes ready for functional and comparative genomics.

## 👤 Author & Contact
**Yasin Kaya**  
Max Planck Institute for Plant Breeding Research  
Email: [ykaya@mpipz.mpg.de](mailto:ykaya@mpipz.mpg.de) | yyasinkkaya@gmail.com

## 📋 Table of Contents
- [Pipeline Overview](#pipeline-overview)  
- [Prerequisites](#prerequisites)  
- [Directory Structure](#directory-structure)  
- [Usage](#usage)  
- [References](#references)  

---

## 🚀 Pipeline Overview

1. **Hybrid Assembly (Verkko2)**
   Generates a diploid-aware assembly by combining ONT and HiFi reads with [Verkko2](https://github.com/marbl/verkko):
   ```bash
   verkko -d <sample> \
     --hifi /path/<sample>.hifi.fastq.gz \
     --nano /path/<sample>.ont.fastq.gz \
     --ref /path/ColCC_reference.fasta \
     --telomere-motif TTTAGGG \
     --grid --local-memory 400 --local-cpus 50
   ```

2. **Hybrid Assembly (hifiasm)**
   Alternative diploid assembly using [hifiasm](https://github.com/chhylp123/hifiasm):
   ```bash
   hifiasm -t 32 -o <output_prefix> \
     --hg-size 140m \
     -ul /path/<sample>.ont.fastq.gz \
     --ul-cut 100000 \
     --dual-scaf --primary --telo-m TTTAGGG \
     /path/<sample>.hifi.fastq.gz
   ```

3. **Nucmer Alignment & QuickMerge**
   Aligns **self** (hifiasm) and **hybrid** (Verkko2) assemblies with [MUMmer4](https://github.com/mummer4/mummer) and merges contigs:
   ```bash
   nucmer -t 10 -p <prefix>.aln self.fasta hybrid.fasta
   delta-filter -r -q -l 10000 -i 0.9 <prefix>.aln.delta > <prefix>.aln.flt.delta
   quickmerge -d <prefix>.aln.flt.delta \
     -q hybrid.fasta -r self.fasta \
     -hco 5.0 -c 1.5 \
     -l <N50_cutoff> -ml 5000 \
     -p <output_prefix> -v
   ```
   **Tip**: Set `-l` to your assembly’s N50 (bp). Smaller values increase merge sensitivity but may introduce mis-joins.

4. **Purge Haplotigs**
   - Generate coverage histograms (`pbcstat`, `calcuts`, `hist_plot.py`).
   - Remove duplicate/haplotig contigs with [Purge_Dups](https://github.com/dfguan/purge_dups).

5. **Scaffolding with RagTag**
   Order and orient contigs against the Col-CC T2T reference (GCA_028009825.2) using [RagTag scaffold](https://github.com/malonge/RagTag).

6. **Gap-Closing**
   Fill remaining gaps with Verkko2-corrected HiFi reads via MaSuRCA’s `close_scaffold_gaps.sh`.

7. **Pilon Polishing**
   Perform three rounds of base-level polishing with Illumina short reads using [Pilon](https://github.com/broadinstitute/pilon).

---

## 📋 Prerequisites

Ensure these tools are in your `PATH`:

- **Hybrid assembly:** `verkko`, `hifiasm`  
- **Alignment & merging:** `nucmer`, `delta-filter`, `quickmerge`  
- **Haplotig purging:** `pbcstat`, `calcuts`, `hist_plot.py`, `split_fa`, `minimap2`, `purge_dups`, `get_seqs`  
- **Scaffolding:** `ragtag.py`, `samtools`, `make_agp.py`  
- **Gap closing:** MaSuRCA (`close_scaffold_gaps.sh`)  
- **Polishing:** `pilon`, `minimap2`, `samtools`

Additional Python packages: `numpy`, `matplotlib`.

---

## 🗂️ Directory Structure

```text
raw_data/                          # Sequencing reads
├── Illumina/<sample>_R1.fastq.gz
├── Illumina/<sample>_R2.fastq.gz
├── HiFi/<sample>.hifi.fastq.gz
└── ONT/<sample>.ont.fastq.gz

genome_assembly/
├── first_assembly/                # Verkko2 & hifiasm outputs
│   ├── verkko2/<sample>/assembly.fasta
│   └── hifiasm/<sample>.asm.p_ctg.fa
├── <sample>/                      # Merged & polished directories
│   ├── 05_quickmerge/
│   ├── 06_scaffolding/
│   ├── 07_gapclosing/
│   └── 08_pilon/
└── Col-CC_reference.fasta          # GCA_028009825.2 T2T reference

AraFusion.sh                       # Main pipeline script
README.md                          # This file
```  

---

## ⚙️ Usage

1. **Customize** file paths, thread counts, and sample IDs at the top of `AraFusion.sh`.  
2. **Run**:
   ```bash
   bash AraFusion.sh
   ```
3. **Inspect**: Check per-sample output folders for logs, plots, and final assemblies.

---

## 📖 References

- Col-CC T2T Reference: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028009825.2/  
- Verkko2: https://github.com/marbl/verkko  
- hifiasm: https://github.com/chhylp123/hifiasm  
- QuickMerge: https://github.com/mahulchak/quickmerge  
- RagTag: https://github.com/malonge/RagTag  
- Pilon: https://github.com/broadinstitute/pilon

---

## License

This project is available under the MIT License.

---

*Good luck!*
