# Conserved Synteny Detection Pipeline

This Python script automates the detection of **conserved synteny** between multiple species using their **GFF3 annotations** and **CDS FASTA files** from [Ensembl](https://useast.ensembl.org/index.html), and a **chromosome mapping file**. It leverages the [JCVI toolkit](https://github.com/tanghaibao/jcvi) and [MCscan-python version](https://github.com/tanghaibao/jcvi/wiki/Mcscan-(python-version)) to perform pairwise orthology detection and generates filtered synteny tables for defined chromosomes.

## Directory Structure

Upon execution, the script organizes intermediate and output files into the following directories:

```
synteny_work/
├── beds/                # BED files derived from GFFs
├── cds/                 # Reformatted CDS FASTA files
├── anchors/             # Raw ortholog anchors + LAST alignment outputs
├── csvs/                # Parsed and merged synteny tables
└── csvs/final_synteny_csvs/  # Final chromosome-filtered synteny results
```

---

## Input Requirements

Place the following files downloaded from [Ensembl](https://useast.ensembl.org/index.html) in the script's parent directory:

### 1. GFF Files

* Place in a folder called: `./gff_files/`
* Format: `.gff3` or `.gff3.gz`
* Filename format: `{species}.gff3(.gz)`

### 2. CDS Files

* Place in a folder called: `./cds_files/`
* Format: `.fa` or `.fa.gz`
* Filename should include the species name to match GFF

### 3. Chromosome Mapping - you make this 

* File: `chromosome_map.json`
* Format:

```json
{
  "species1": "chrX",
  "species2": "ChrY",
  ...
}
```

This dictionary is used to filter orthologous pairs that fall within specified chromosome pairs. An example file is available for you. 

---

## Dependencies

Ensure the following are installed and available in your environment:

* Python ≥ 3.6
* [JCVI](https://github.com/tanghaibao/jcvi) (`pip install jcvi`)
* `pandas`
* `gzip`, `subprocess`, `shutil`, `os` — standard Python libraries

---

## Script Overview

### 1. **Preparation**

* Create working directories.
* Load chromosome map.
* Identify matching GFF and CDS files by species name.

### 2. **Preprocessing**

* Convert GFF to BED format using `jcvi.formats.gff`.
* Reformat CDS using `jcvi.formats.fasta`.
* Normalize gene IDs (remove prefixes and version numbers).
* Check overlap between BED and CDS IDs for consistency.

### 3. **Ortholog Search**

* For each species pair:

  * Run `jcvi.compara.catalog ortholog` to detect orthologs.
  * Store anchor, LAST, and PDF outputs in pair-specific folders.

### 4. **Postprocessing**

* Parse `.anchors` files and merge with BED data.
* Optionally replace Ensembl IDs with gene names.
* Output to `synteny_*.csv` files.

### 5. **Filtering**

* Use the chromosome map to filter results for chromosome-specific synteny.
* Save results in `csvs/final_synteny_csvs`.

---

## Usage

```bash
python synteny_pipeline.py
```

The script will automatically:

* Detect all species pairs
* Run synteny analysis
* Save final `.csv` outputs of conserved synteny by chromosome

---

## Output

Filtered synteny tables (one per species pair) are saved in:

```
synteny_work/csvs/final_synteny_csvs/
```

Each file contains:

* Gene pairs
* Alignment score
* Chromosomal positions in both species

---

## Example Output (TSV)

```
gene1   gene2   score   chrom1   start1   end1   chrom2   start2   end2
BRCA1   BRCA1B  972     chr1     123456   125000 chr2     112000   113400
...
```

---

## Notes and Troubleshooting

* If `.anchors` files are missing, the species may have poor LAST alignment or naming mismatches.
* Ensure species names are consistent across GFF, CDS, and chromosome map.
* You can manually inspect the `anchors/` directory for intermediate logs and alignment files.

---
## Plotting
* To create a large circular plot of conserved synteny of your chromosomal regions see: [synteny_rplot](https://github.com/ian-bda/synteny_rplot)

---
Script developed by **Ian Birchler De Allende**
Uses [JCVI Comparative Genomics toolkit](https://github.com/tanghaibao/jcvi) and [MCscan-python version](https://github.com/tanghaibao/jcvi/wiki/Mcscan-(python-version))

---
