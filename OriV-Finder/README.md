# OriV-Finder Batch Runner

This folder contains the OriV-Finder Docker image archive and a batch runner for predicting oriV/Rep features from plasmid FASTA files.

## Folder Layout

```text
OriV-Finder/
  orivfinder-ready.tar.gz
  run_batch_orivfinder.ps1
  data/
    input/      put .fasta, .fa, or .fna files here
    output/     one output folder per input FASTA
```

## 1. Load Docker Image

If the image has not been loaded yet:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_batch_orivfinder.ps1 -LoadImage -DryRun
```

Or load manually after extracting the archive:

```powershell
docker load -i orivfinder-ready.tar
```

Check the image:

```powershell
docker images | findstr orivfinder-ready
```

## 2. Prepare Input

Put plasmid nucleotide FASTA files into:

```text
OriV-Finder/data/input/
```

Supported extensions:

```text
.fasta
.fa
.fna
```

## 3. Run Batch Prediction

From the `OriV-Finder` folder:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_batch_orivfinder.ps1
```

The script runs one Docker job per FASTA file:

```text
python oriVfinder.py --fasta /app/data/input/<file> --output_dir /app/data/output/<file-stem>
```

## 4. Useful Options

Preview commands without running Docker:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_batch_orivfinder.ps1 -DryRun
```

Use a custom input/output folder:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_batch_orivfinder.ps1 ^
  -InputDir data/input ^
  -OutputDir data/output
```

Overwrite completed outputs:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_batch_orivfinder.ps1 -Overwrite
```

Use a different image tag:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_batch_orivfinder.ps1 ^
  -ImageName orivfinder-ready ^
  -ImageTag latest
```

## 5. Output

Each input FASTA gets its own output folder:

```text
data/output/<input-file-stem>/
```

A batch summary is written to:

```text
data/output/batch_summary.tsv
```

Completed inputs get a marker file:

```text
data/output/<input-file-stem>/.done
```

If `.done` exists, the script skips that input unless `-Overwrite` is used.
