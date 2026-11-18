# 🧬 Genome Annotation Parser

A Streamlit web app for parsing, exploring, and visualizing genome annotation
files in **GFF3** or **GTF** format.

## Features

- **Format auto-detection** — handles GFF3 (`key=value`) and GTF (`key "value"`)
  attribute styles, plain text or gzip-compressed (`.gz`).
- **Overview dashboard** — feature type counts, features per chromosome,
  strand distribution, source breakdown.
- **Feature explorer** — filter by type, sequence/chromosome, source, length
  range, or free-text search across ID/Name/attributes; inspect any single
  feature's full attribute set.
- **Visualizations** — length-distribution histograms, a genomic feature map
  (features plotted along a chosen chromosome), and feature-density plots.
- **Export** — download the filtered feature set as CSV or GFF3, or export
  the full parsed dataset.

A small bundled sample GFF3 file (`sample_data/sample.gff3`) lets you try the
app immediately without uploading your own data.


## Input format notes

The app expects the standard 9-column tab-separated GFF3/GTF layout:

```
seqid  source  type  start  end  score  strand  phase  attributes
```

- Lines starting with `#` are treated as comments/headers and skipped.
- A `##FASTA` directive (or embedded sequence block) is detected and
  everything after it is ignored, so you can safely upload GFF3 files with
  an appended FASTA section.
- GFF3-style (`ID=x;Parent=y`) and GTF-style (`gene_id "x"; transcript_id "y";`)
  attribute strings are both parsed into a dictionary you can inspect per
  feature.
