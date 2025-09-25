# Cell Ranger Setup and Run Guide

This document summarizes the useful commands from modifying the `PATH` variable to renaming FASTQ files and running `cellranger count`.

---

## 1. Set PATH (for current session)

```bash
export PATH="/home/qian@vitra.bio/cellranger-8.0.1:$PATH"

# Verify installation
which cellranger
cellranger --version
```

---

## 2. Make PATH permanent

```bash
echo 'export PATH="/home/qian@vitra.bio/cellranger-8.0.1:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

**Optional:** If you prefer moving the folder to a path without `@` (safer for scripts):

```bash
mv "/home/qian@vitra.bio/cellranger-8.0.1" "/home/qian/cellranger-8.0.1"
export PATH="/home/qian/cellranger-8.0.1:$PATH"
```

---

## 3. Rename FASTQ files to Cell Ranger naming convention

Inside the sample folder `SRR33132543`:

```bash
cd "/home/qian@vitra.bio/ava/data/huang_data_mouse_whole_ovary/single_cell/SRR33132543"

mv SRR33132543_1.fastq.gz SRR33132543_S1_L001_I1_001.fastq.gz
mv SRR33132543_2.fastq.gz SRR33132543_S1_L001_I2_001.fastq.gz
mv SRR33132543_3.fastq.gz SRR33132543_S1_L001_R1_001.fastq.gz
mv SRR33132543_4.fastq.gz SRR33132543_S1_L001_R2_001.fastq.gz

ls -lh
```

---

## 4. Run `cellranger count`

```bash
cellranger count   --id="SRR33132543_out"   --transcriptome="/home/qian@vitra.bio/refdata-gex-GRCm39-2024-A"   --fastqs="/home/qian@vitra.bio/ava/data/huang_data_mouse_whole_ovary/single_cell/SRR33132543"   --sample="SRR33132543"   --create-bam=false
```

The output will be in `SRR33132543_out/`.  
Check logs with:

```bash
tail -f SRR33132543_out/_log
```

---

## 5. Batch script for multiple SRR samples

Save this as `run_cellranger_batch.sh` and run with `bash run_cellranger_batch.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

TRANSCRIPTOME="/home/qian@vitra.bio/refdata-gex-GRCm39-2024-A"
BASE_DIR="/home/qian@vitra.bio/ava/data/huang_data_mouse_whole_ovary/single_cell"

for d in "${BASE_DIR}"/SRR*; do
  [ -d "$d" ] || continue
  sample=$(basename "$d")

  echo ">>> Processing ${sample}"

  if ls "${d}/${sample}_S1_L001_R1_001.fastq.gz" >/dev/null 2>&1; then
    echo "    - FASTQ already renamed. Skip renaming."
  else
    mv "${d}/${sample}_1.fastq.gz" "${d}/${sample}_S1_L001_I1_001.fastq.gz"
    mv "${d}/${sample}_2.fastq.gz" "${d}/${sample}_S1_L001_I2_001.fastq.gz"
    mv "${d}/${sample}_3.fastq.gz" "${d}/${sample}_S1_L001_R1_001.fastq.gz"
    mv "${d}/${sample}_4.fastq.gz" "${d}/${sample}_S1_L001_R2_001.fastq.gz"
    echo "    - Renamed FASTQs to Cell Ranger convention."
  fi

  outdir="${BASE_DIR}/${sample}_out"
  if [ -d "$outdir" ]; then
    echo "    - Output exists (${outdir}), skip running."
    continue
  fi

  echo "    - Running cellranger count..."
  cellranger count     --id="${sample}_out"     --transcriptome="${TRANSCRIPTOME}"     --fastqs="${d}"     --sample="${sample}"     --create-bam=false

  echo ">>> Done ${sample}"
done
```
