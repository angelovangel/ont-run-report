## TGS Run Report

An interactive Shiny web application that generates an HTML report for Oxford Nanopore Technologies (ONT) or PacBio HiFi sequencing runs.

---

### How It Works

Select the sequencer type, upload the required CSV statistics files, configure a few options, and click **Generate Report**. The app calls the appropriate Python script in the background and when finished, an **Open / Download HTML Report** button appears.

---

### ONT mode

Select the two MinKNOW CSV files for every flow cell in a single upload; they
are paired automatically by the run ID in their file names, and each flow cell
gets a panel where you can change its label (default: the flow cell ID). Files
that can't be paired are listed and ignored. Sequencing summary files are
optional; add one in every flow cell's panel to include the interactive
read-length Nx comparison:

| File | Description |
| :--- | :--- |
| `pore_activity_*.csv` | Per-pore-state time series (number of pores in each state over the run). |
| `throughput_*.csv` | Cumulative yield / read-count time series. |
| `sequencing_summary.txt` / TSV | Optional local file for the Nx plot (`sequence_length_template` column); sampled in the browser.

These files are found inside the run directory produced by MinKNOW, typically under `<experiment>/<run_id>/`. They are named `pore_activity_<run_id>.csv`, `throughput_<run_id>.csv` and `sequencing_summary_<run_id>.txt`.

---

#### Command line

The report script also runs on its own (Python 3, no extra packages). Give it
MinKNOW run folders and/or files in any order; files are matched by the run ID
in their names and folders are searched recursively:

```bash
bin/generate-ont-report.py /data/run1 /data/run2 \
    -l PBC87221=Liver PBE26402=Blood -t "Week 42" -o week42.html
```

With no inputs it searches the current folder. Labels default to the flow cell
ID; `-l` takes `FLOWCELL=label` (or `RUNID=label`), or plain labels in the
order the runs are listed. The Nx plot is added when every run has a
`sequencing_summary_*.txt`. See `bin/generate-ont-report.py -h` for all options.

---

### PacBio HiFi mode

Upload one or more `*run-qc-export*.csv` files exported from **SMRT Link** (one file per run / instrument). Multiple files are merged into a single report, which includes:

- HiFi yield per sample, with each SMRT cell as a segment of the sample's bar;
- a cell table grouped by run, with run subtotals and a total that follow the
  filters; read lengths (polymerase RL, N50, insert size, mean length) and
  loading (P0 / P1 / P2) are shown as small inline charts, and the table
  downloads as CSV;
- any metric per cell, metric-vs-metric correlation, and — when a SMRT cell
  was used more than once — a metric across its acquisitions;
- optionally, a **Demultiplexing** section: add SMRT Link
  `*.barcodes.summary.csv` files to the upload (or pass them with `--demux`)
  to get reads, yield, read length, QV and barcode quality per barcode, with
  each barcode's share of the cell's reads and the CV across barcodes. Each
  file is matched to its acquisition by the instrument, movie date and sample
  name in its file name (`analysis-<sample>-<job>-m<instrument>_<YYMMDD>_<HHMMSS>_s<N>…`),
  then by the HiFi yield closest to the barcodes' total.

---

> Note: the provided Dockerfile is for deployment with [tgs-lab-apps](https://github.com/angelovangel/tgs-lab-apps).
