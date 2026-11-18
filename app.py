"""
Genome Annotation Parser

"""

import io
import re
import gzip
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# --------------------------------------------------------------------------
# Page config
# --------------------------------------------------------------------------
st.set_page_config(
    page_title="Genome Annotation Parser",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

GFF_COLUMNS = [
    "seqid", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes",
]


# --------------------------------------------------------------------------
# Parsing helpers
# --------------------------------------------------------------------------
def _read_raw_text(uploaded_file) -> str:
    """Read an uploaded file (optionally gzip-compressed) into a text string."""
    raw_bytes = uploaded_file.read()
    if uploaded_file.name.endswith(".gz"):
        raw_bytes = gzip.decompress(raw_bytes)
    return raw_bytes.decode("utf-8", errors="replace")


def _detect_format(sample_lines) -> str:
    """Guess whether the attribute column looks like GFF3 (key=value) or GTF (key "value")."""
    for line in sample_lines:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 9:
            attr = parts[8]
            if re.search(r'\w+\s+"', attr):
                return "GTF"
            if re.search(r"\w+=", attr):
                return "GFF3"
    return "GFF3"


def _parse_attributes(attr_string: str, fmt: str) -> dict:
    """Turn the raw attribute column into a dict of key -> value."""
    result = {}
    attr_string = attr_string.strip().rstrip(";")
    if not attr_string:
        return result

    if fmt == "GFF3":
        for field in attr_string.split(";"):
            field = field.strip()
            if not field or "=" not in field:
                continue
            key, _, value = field.partition("=")
            result[key.strip()] = value.strip()
    else:  # GTF
        for field in attr_string.split(";"):
            field = field.strip()
            if not field:
                continue
            m = re.match(r'(\S+)\s+"?([^"]*)"?', field)
            if m:
                key, value = m.group(1), m.group(2)
                result[key.strip()] = value.strip()
    return result


@st.cache_data(show_spinner=False)
def parse_annotation_file(file_bytes: bytes, file_name: str):
    """Parse a GFF3/GTF file (given as raw bytes) into a tidy DataFrame."""
    if file_name.endswith(".gz"):
        text = gzip.decompress(file_bytes).decode("utf-8", errors="replace")
    else:
        text = file_bytes.decode("utf-8", errors="replace")

    lines = text.splitlines()
    fmt = _detect_format(lines[:500])

    records = []
    skipped = 0
    fasta_started = False
    for line in lines:
        if not line or line.startswith("##FASTA"):
            fasta_started = True
        if fasta_started:
            break
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            skipped += 1
            continue
        seqid, source, ftype, start, end, score, strand, phase, attrs = parts[:9]
        try:
            start_i = int(start)
            end_i = int(end)
        except ValueError:
            skipped += 1
            continue

        attr_dict = _parse_attributes(attrs, fmt)
        record_id = (
            attr_dict.get("ID")
            or attr_dict.get("gene_id")
            or attr_dict.get("transcript_id")
            or attr_dict.get("Name")
            or ""
        )
        parent = attr_dict.get("Parent") or attr_dict.get("gene_id") or ""
        name = attr_dict.get("Name") or attr_dict.get("gene_name") or record_id

        records.append({
            "seqid": seqid,
            "source": source,
            "type": ftype,
            "start": start_i,
            "end": end_i,
            "length": end_i - start_i + 1,
            "score": score,
            "strand": strand,
            "phase": phase,
            "id": record_id,
            "name": name,
            "parent": parent,
            "attributes": attrs,
            "attr_dict": attr_dict,
        })

    df = pd.DataFrame.from_records(records)
    if not df.empty:
        df["score_num"] = pd.to_numeric(df["score"], errors="coerce")
    return df, fmt, skipped


def df_to_gff3(df: pd.DataFrame) -> str:
    """Serialize a (filtered) DataFrame back into GFF3 text."""
    lines = ["##gff-version 3"]
    for _, row in df.iterrows():
        lines.append("\t".join([
            str(row["seqid"]), str(row["source"]), str(row["type"]),
            str(row["start"]), str(row["end"]), str(row["score"]),
            str(row["strand"]), str(row["phase"]), str(row["attributes"]),
        ]))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Sidebar - file input
# --------------------------------------------------------------------------
st.sidebar.title("🧬 Genome Annotation Parser")
st.sidebar.markdown(
    "Upload a **GFF3** or **GTF** annotation file (`.gff`, `.gff3`, `.gtf`, "
    "optionally gzip-compressed as `.gz`)."
)

uploaded_file = st.sidebar.file_uploader(
    "Annotation file",
    type=["gff", "gff3", "gtf", "gz"],
    help="Standard 9-column GFF3/GTF format.",
)

use_sample = st.sidebar.checkbox("Use bundled sample data instead", value=False)

if use_sample:
    try:
        with open("sample_data/sample.gff3", "rb") as f:
            file_bytes = f.read()
        file_name = "sample.gff3"
    except FileNotFoundError:
        st.sidebar.error("Sample file not found in sample_data/.")
        st.stop()
elif uploaded_file is not None:
    file_bytes = uploaded_file.read()
    file_name = uploaded_file.name
else:
    st.title("🧬 Genome Annotation Parser")
    st.info(
        "👈 Upload a GFF3/GTF file from the sidebar, or check "
        "**'Use bundled sample data'** to try the app right away."
    )
    st.markdown(
        """
        ### What this app does
        - Parses **GFF3** and **GTF** genome annotation files
        - Summarizes feature types, sources, and sequence regions
        - Lets you filter and search features interactively
        - Visualizes feature distribution, lengths, and genomic position
        - Exports filtered results as **CSV** or **GFF3**
        """
    )
    st.stop()

with st.spinner("Parsing annotation file..."):
    df, detected_format, skipped_lines = parse_annotation_file(file_bytes, file_name)

if df.empty:
    st.error(
        "No valid annotation records were found. Please check that the file "
        "is a properly formatted GFF3/GTF file."
    )
    st.stop()

st.sidebar.success(f"Detected format: **{detected_format}**")
st.sidebar.write(f"**{len(df):,}** features parsed")
if skipped_lines:
    st.sidebar.warning(f"{skipped_lines} malformed line(s) skipped")

# --------------------------------------------------------------------------
# Sidebar - filters
# --------------------------------------------------------------------------
st.sidebar.markdown("---")
st.sidebar.subheader("Filters")

all_types = sorted(df["type"].unique().tolist())
selected_types = st.sidebar.multiselect("Feature type", all_types, default=[])

all_seqids = sorted(df["seqid"].unique().tolist())
selected_seqids = st.sidebar.multiselect("Sequence / chromosome", all_seqids, default=[])

all_sources = sorted(df["source"].unique().tolist())
selected_sources = st.sidebar.multiselect("Source", all_sources, default=[])

search_term = st.sidebar.text_input("Search ID / Name / attributes contains…")

min_len, max_len = int(df["length"].min()), int(df["length"].max())
if min_len < max_len:
    length_range = st.sidebar.slider(
        "Feature length range", min_value=min_len, max_value=max_len,
        value=(min_len, max_len),
    )
else:
    length_range = (min_len, max_len)

filtered = df.copy()
if selected_types:
    filtered = filtered[filtered["type"].isin(selected_types)]
if selected_seqids:
    filtered = filtered[filtered["seqid"].isin(selected_seqids)]
if selected_sources:
    filtered = filtered[filtered["source"].isin(selected_sources)]
if search_term:
    mask = (
        filtered["id"].str.contains(search_term, case=False, na=False)
        | filtered["name"].str.contains(search_term, case=False, na=False)
        | filtered["attributes"].str.contains(search_term, case=False, na=False)
    )
    filtered = filtered[mask]
filtered = filtered[
    (filtered["length"] >= length_range[0]) & (filtered["length"] <= length_range[1])
]

# --------------------------------------------------------------------------
# Main - tabs
# --------------------------------------------------------------------------
st.title("🧬 Genome Annotation Parser")
st.caption(f"File: `{file_name}` · Format: {detected_format}")

tab_overview, tab_explorer, tab_viz, tab_export = st.tabs(
    ["📊 Overview", "🔍 Feature Explorer", "📈 Visualizations", "⬇️ Export"]
)

# ---- Overview tab ---------------------------------------------------------
with tab_overview:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total features", f"{len(df):,}")
    c2.metric("Feature types", f"{df['type'].nunique():,}")
    c3.metric("Sequences / chromosomes", f"{df['seqid'].nunique():,}")
    c4.metric("Sources", f"{df['source'].nunique():,}")

    st.subheader("Feature type counts")
    type_counts = df["type"].value_counts().reset_index()
    type_counts.columns = ["type", "count"]
    fig = px.bar(type_counts, x="type", y="count", color="type",
                 title="Number of features by type")
    fig.update_layout(showlegend=False, xaxis_title="", yaxis_title="Count")
    st.plotly_chart(fig, use_container_width=True)

    col_a, col_b = st.columns(2)
    with col_a:
        st.subheader("Features per sequence")
        seq_counts = df["seqid"].value_counts().reset_index()
        seq_counts.columns = ["seqid", "count"]
        fig2 = px.bar(seq_counts.head(30), x="seqid", y="count",
                      title="Top sequences by feature count")
        fig2.update_layout(xaxis_title="", yaxis_title="Count")
        st.plotly_chart(fig2, use_container_width=True)

    with col_b:
        st.subheader("Strand distribution")
        strand_counts = df["strand"].value_counts().reset_index()
        strand_counts.columns = ["strand", "count"]
        fig3 = px.pie(strand_counts, names="strand", values="count",
                      title="Features by strand")
        st.plotly_chart(fig3, use_container_width=True)

    st.subheader("Sources")
    st.dataframe(df["source"].value_counts().rename("count"), use_container_width=True)

# ---- Feature Explorer tab -------------------------------------------------
with tab_explorer:
    st.write(f"Showing **{len(filtered):,}** of **{len(df):,}** features "
             "after filters (use the sidebar to adjust).")

    display_cols = ["seqid", "source", "type", "start", "end", "length",
                     "strand", "score", "id", "name", "parent"]
    st.dataframe(filtered[display_cols], use_container_width=True, height=420)

    st.markdown("#### Inspect a single feature")
    if len(filtered) > 0:
        idx = st.number_input(
            "Row index (from the table above, 0-based)", min_value=0,
            max_value=max(len(filtered) - 1, 0), value=0, step=1,
        )
        row = filtered.iloc[int(idx)]
        colx, coly = st.columns([1, 1])
        with colx:
            st.json({
                "seqid": row["seqid"], "source": row["source"], "type": row["type"],
                "start": int(row["start"]), "end": int(row["end"]),
                "length": int(row["length"]), "strand": row["strand"],
                "score": row["score"], "phase": row["phase"],
            })
        with coly:
            st.json(row["attr_dict"])

# ---- Visualizations tab ---------------------------------------------------
with tab_viz:
    st.subheader("Feature length distribution")
    type_for_len = st.selectbox(
        "Restrict to feature type (optional)",
        ["All"] + all_types, index=0, key="len_type",
    )
    len_df = df if type_for_len == "All" else df[df["type"] == type_for_len]
    fig4 = px.histogram(len_df, x="length", nbins=60,
                         title=f"Length distribution ({type_for_len})")
    fig4.update_layout(xaxis_title="Length (bp)", yaxis_title="Count")
    st.plotly_chart(fig4, use_container_width=True)

    st.subheader("Genomic feature map")
    seq_choice = st.selectbox("Choose a sequence / chromosome", all_seqids, key="seq_map")
    seq_df = df[df["seqid"] == seq_choice].copy()

    type_choices = st.multiselect(
        "Feature types to show on the map", all_types,
        default=all_types[: min(5, len(all_types))],
    )
    if type_choices:
        seq_df = seq_df[seq_df["type"].isin(type_choices)]

    if not seq_df.empty:
        color_map = {t: px.colors.qualitative.Set2[i % 8]
                     for i, t in enumerate(sorted(seq_df["type"].unique()))}
        fig5 = go.Figure()
        y_positions = {t: i for i, t in enumerate(sorted(seq_df["type"].unique()))}
        for t in sorted(seq_df["type"].unique()):
            sub = seq_df[seq_df["type"] == t]
            fig5.add_trace(go.Bar(
                x=sub["length"], y=[t] * len(sub), base=sub["start"],
                orientation="h", name=t, marker_color=color_map[t],
                hovertext=sub.apply(
                    lambda r: f"{r['name']} ({r['start']}-{r['end']}, {r['strand']})",
                    axis=1,
                ),
                hoverinfo="text",
            ))
        fig5.update_layout(
            title=f"Feature map for {seq_choice}",
            xaxis_title="Position (bp)", yaxis_title="Feature type",
            barmode="overlay", height=120 + 40 * len(y_positions),
        )
        st.plotly_chart(fig5, use_container_width=True)
    else:
        st.info("No features match the selected types on this sequence.")

    st.subheader("Feature density along sequence")
    if not seq_df.empty:
        bins = st.slider("Number of bins", 10, 200, 50)
        seq_len = df[df["seqid"] == seq_choice]["end"].max()
        bin_edges = np.linspace(0, seq_len, bins + 1)
        seq_df["bin"] = pd.cut(seq_df["start"], bins=bin_edges)
        density = seq_df.groupby("bin", observed=True).size().reset_index(name="count")
        density["bin_start"] = density["bin"].apply(lambda b: b.left).astype(float)
        fig6 = px.area(density, x="bin_start", y="count",
                        title=f"Feature density across {seq_choice}")
        fig6.update_layout(xaxis_title="Position (bp)", yaxis_title="Feature count")
        st.plotly_chart(fig6, use_container_width=True)

# ---- Export tab -------------------------------------------------------------
with tab_export:
    st.write("Export the **currently filtered** feature set (see Feature Explorer filters).")
    st.write(f"{len(filtered):,} features will be exported.")

    csv_bytes = filtered[display_cols].to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download filtered features as CSV",
        data=csv_bytes,
        file_name="filtered_annotations.csv",
        mime="text/csv",
    )

    gff3_text = df_to_gff3(filtered)
    st.download_button(
        "Download filtered features as GFF3",
        data=gff3_text.encode("utf-8"),
        file_name="filtered_annotations.gff3",
        mime="text/plain",
    )

    st.markdown("---")
    st.write("Export full parsed dataset (no filters applied).")
    full_csv = df[display_cols].to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download full dataset as CSV",
        data=full_csv,
        file_name="all_annotations.csv",
        mime="text/csv",
    )

st.sidebar.markdown("---")
st.sidebar.caption("Built with Streamlit · pandas · plotly")
