## ONT Run Report

An interactive Shiny web application that generates an HTML quality-control report for one or more Oxford Nanopore Technologies (ONT) sequencing flow cells.

---

### How It Works

Upload the CSV statistics files exported by MinKNOW for each flow cell, configure a few options, and click **Generate Report**. The app calls `generate-ont-report.py` in the background and streams live terminal output. When finished, a **Download HTML Report** button appears.

---

### Input Files

For each flow cell you need two CSV files (exported from the MinKNOW run directory):

| File | Description |
| :--- | :--- |
| `pore_activity_*.csv` | Per-pore-state time series (number of pores in each state over the run). |
| `throughput_*.csv` | Cumulative yield / read-count time series. |

Both files are found inside the run directory produced by MinKNOW, typically under `<experiment>/<run_id>/`. They are named `pore_activity_<run_id>.csv` and `throughput_<run_id>.csv`.

---

### Options

| Option | Default | Description |
| :--- | :--- | :--- |
| **Number of flow cells** | 2 | How many flow cells to compare in the report (1-8). |
| **Label** | Run N | Human-readable name for each flow cell. |
| **Report title** | ONT Sequencing Report | Title shown at the top of the HTML report. |
| **Active-pores sampling (min)** | 5 | Down-sampling interval for the Active-Pores chart. |

---

### UI Features

* **Dynamic file inputs**: Set the number of flow cells; the app renders the corresponding file-input groups automatically.
* **Command preview**: Before running, the sidebar terminal shows the exact Python command that will be executed.
* **Live log streaming**: Once started, the terminal streams stdout/stderr in real time.
* **Download button**: Appears after a successful run; click to save the self-contained HTML report.
* **Kill process**: Safely terminates a running job at any time.
* **Reset**: Clears all inputs and temporary files.

---

### Getting Started

#### 1. Interactive Web App (Docker)

```bash
docker compose pull
docker compose up -d ont-run-report
```

Then open `http://<server>:3810` in your browser.

#### 2. Command Line

```bash
python generate-ont-report.py \
  --pore-activity pore_activity_run1.csv pore_activity_run2.csv \
  --throughput    throughput_run1.csv     throughput_run2.csv \
  --labels        "Run 1" "Run 2" \
  --title         "My Experiment" \
  --out           report.html
```
