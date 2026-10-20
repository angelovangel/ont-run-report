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

### PacBio HiFi mode

Upload one or more `*run-qc-export*.csv` files exported from **SMRT Link** (one file per run / instrument). Multiple files are merged into a single report, which includes:

- HiFi yield per sample, with each SMRT cell as a segment of the sample's bar;
- a cell table grouped by run, with run subtotals and a total that follow the
  filters; read lengths (polymerase RL, N50, insert size, mean length) and
  loading (P0 / P1 / P2) are shown as small inline charts, and the table
  downloads as CSV;
- any metric per cell, metric-vs-metric correlation, and — when a SMRT cell
  was used more than once — a metric across its acquisitions.

---

> Note: the provided Dockerfile is for deployment with [tgs-lab-apps](https://github.com/angelovangel/tgs-lab-apps).
