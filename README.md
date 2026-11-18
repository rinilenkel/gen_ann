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

## Project structure

```
genome-annotation-parser/
├── app.py                     # main Streamlit app
├── requirements.txt           # Python dependencies
├── sample_data/
│   └── sample.gff3            # demo annotation file
├── .streamlit/
│   └── config.toml            # theme + upload size config
└── README.md
```

## Running locally

Requires Python 3.9+.

```bash
# 1. Clone your repo (after you've pushed this project to GitHub)
git clone https://github.com/<your-username>/genome-annotation-parser.git
cd genome-annotation-parser

# 2. (Recommended) create a virtual environment
python -m venv venv
source venv/bin/activate      # Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Run the app
streamlit run app.py
```

The app will open at `http://localhost:8501`. Check **"Use bundled sample
data"** in the sidebar to explore it right away, or upload your own
`.gff`, `.gff3`, `.gtf`, or gzip-compressed variant.

## Deploying to GitHub + Streamlit Community Cloud

1. **Create a GitHub repository** and push this project:

   ```bash
   git init
   git add .
   git commit -m "Initial commit: genome annotation parser app"
   git branch -M main
   git remote add origin https://github.com/<your-username>/genome-annotation-parser.git
   git push -u origin main
   ```

2. **Deploy on Streamlit Community Cloud**:
   - Go to [share.streamlit.io](https://share.streamlit.io) and sign in with
     GitHub.
   - Click **"New app"**.
   - Select your repository, branch (`main`), and set the main file path to
     `app.py`.
   - Click **"Deploy"**.

   Streamlit Cloud will automatically install everything listed in
   `requirements.txt` and launch the app. Any future `git push` to `main`
   will trigger an automatic redeploy.

3. **Optional settings**:
   - If you expect large annotation files, you can raise
     `server.maxUploadSize` in `.streamlit/config.toml` (already set to
     500 MB) — Streamlit Cloud enforces its own platform limits on top of
     this, so very large genomes may need to be pre-filtered locally.

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

## Extending the app

Some ideas if you want to build on this:

- Add BLAST/InterPro annotation cross-referencing
- Support multi-file comparison (e.g., two annotation versions/assemblies)
- Add codon/GC-content analysis using a companion FASTA file
- Swap the custom parser for `gffutils` or `BioPython.SeqIO` if you need
  full GFF3 spec compliance (hierarchical feature graphs, validation, etc.)

## License

Use and adapt freely for your own research or teaching purposes.
