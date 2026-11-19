# 🧬 Genome Annotation Parser

An interactive web application built with **Python**, **Biopython**, and **Streamlit** to parse, filter, visualize, and export genomic feature annotations from **GenBank (`.gb` / `.gbk`)** files. 

---

## 📋 Overview

Raw genome annotation files contain rich biological metadata—such as gene boundaries, coding sequences (CDS), non-coding RNAs, and functional protein descriptions—that can be difficult to inspect manually. This tool converts raw GenBank records into a structured tabular format and renders an interactive genomic coordinate map directly in the browser.

### ✨ Key Features
* **Automated GenBank Parsing:** Uses `Bio.SeqIO` to extract sequence-level metrics and feature-level qualifiers in memory without writing temporary files to disk.
* **Genome Summary Dashboard:** Calculates and displays Accession ID, total genome size (bp), overall GC content percentage (`Bio.SeqUtils.gc_fraction`), and total annotated feature count.
* **Dynamic Filtering & Search:** Filter features by biological category (`gene`, `CDS`, `tRNA`, `rRNA`, etc.) and perform real-time keyword searches across gene symbols, locus tags, and product descriptions.
* **Interactive Feature Track Map:** Renders a zoomable, horizontal coordinate plot using **Plotly Express** (`px.timeline`) to visualize feature positions and strand orientation (`+` / `-`) along the replicon.
* **One-Click CSV Export:** Download filtered annotation tables directly as `.csv` files for downstream bioinformatics analysis.

---

## 🛠️ Tech Stack & Dependencies

| Library | Purpose |
| :--- | :--- |
| **Python 3.10+** | Core programming language |
| **Biopython (`Bio`)** | GenBank file parsing (`SeqIO`) and GC content computation (`gc_fraction`) |
| **Streamlit** | Web application framework and interactive UI components |
| **Pandas** | Tabular feature structuring, boolean filtering, and CSV serialization |
| **Plotly Express** | Interactive genomic coordinate visualization |

---

## 📂 Repository Structure

```text
genome-annotation-parser/
├── app.py               # Main Streamlit web application and Biopython parser
├── requirements.txt     # Python dependencies for local setup & cloud deployment
├── sample_data/         # Optional folder containing sample .gbk files for testing
└── README.md            # Project documentation
