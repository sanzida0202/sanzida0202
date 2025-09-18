# Hi, I'm @sanzida! 👋

I'm a molecular diagnostics enthusiast who loves translating bioinformatics analyses into actionable assay insights. Below is a featured project that now powers my GitHub profile.

## Microsatellite Instability vs. Immunohistochemistry Concordance

**Microsatellite Instability vs Immunohistochemistry Concordance** is a reproducible analysis that evaluates how well next-generation sequencing (NGS) derived MSI calls align with mismatch repair (MMR) protein staining by immunohistochemistry (IHC) in colorectal cancer samples.

### Project Overview
- **Goal:** Quantify concordance and discordance between MSI status (from NGS) and IHC-based MMR interpretation to support assay validation.
- **Implementation:** Python with `pandas` for data wrangling, `statsmodels` for logistic regression, and `matplotlib` for publication-ready visualisations.
- **Highlights:**
  - Automated data ingestion for MSI calls and IHC marker results.
  - Calculation of concordance rates, discordance taxonomy, and contingency tables.
  - Threshold optimisation using Youden's J statistic to propose an updated MSI score cut-off.
  - Logistic regression modelling that translates continuous MSI scores into the probability of IHC-defined dMMR.
  - Exportable figures and CSV/JSON summaries for inclusion in validation documentation.

### Repository Structure
```
├── data/
│   ├── ihc_results.csv          # Example IHC marker staining calls (retained/lost)
│   └── msi_calls.csv            # Example MSI scores and historical MSI calls
├── msi_ihc_concordance/
│   ├── __init__.py
│   └── analysis.py              # Core analysis helpers (data prep, stats, plotting)
├── outputs/                     # Target location for generated reports/figures
├── scripts/
│   └── run_msi_concordance_analysis.py  # CLI entry point for the workflow
├── requirements.txt
└── README.md
```

### Running the Analysis Locally
1. **Create an isolated environment (recommended):**
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\\Scripts\\activate
   ```
2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```
3. **Execute the analysis workflow:**
   ```bash
   python scripts/run_msi_concordance_analysis.py \
       --msi-data data/msi_calls.csv \
       --ihc-data data/ihc_results.csv \
       --output-dir outputs
   ```
   The script prints a JSON summary to the console and writes detailed outputs (merged tables, discordant case list, summary JSON, and a figure) to the `outputs/` directory.

### Sample Output
Using the bundled example data, the workflow produces the following high-level summary:
```json
{
  "total_cases": 12,
  "concordant_cases": 9,
  "discordant_cases": 3,
  "concordance_rate": 0.75,
  "discordance_rate": 0.25,
  "contingency_table": {
    "MSI-H": {"dMMR": 5, "pMMR": 1},
    "MSS": {"pMMR": 4, "dMMR": 2}
  },
  "discordant_case_ids": ["CRC004", "CRC009", "CRC010"],
  "optimal_threshold": {
    "threshold": 21.6,
    "youden_index": 0.71,
    "sensitivity": 0.71,
    "specificity": 1.0
  },
  "logistic_regression": {
    "intercept": -3.25,
    "msi_score_coef": 0.17,
    "p_values": {
      "intercept": 0.13,
      "msi_score": 0.10
    }
  }
}
```
*Values are rounded for readability; running the code will output the full precision metrics as well as CSV and figure artifacts.*

### Generated Visualisations and Reports
- `outputs/msi_score_by_ihc_status.png`: Box/strip plot illustrating MSI score distributions stratified by IHC status with the proposed threshold overlay.
- `outputs/merged_results.csv`: Combined dataset containing MSI calls, IHC interpretations, concordance flag, and logistic regression probabilities.
- `outputs/discordant_cases.csv`: Focused table of cases requiring orthogonal review.
- `outputs/concordance_summary.json`: Machine-readable summary of all key statistics for traceability in validation reports.

### Extending the Project
- Replace the example CSV files with your laboratory's datasets to reproduce the analysis on new cohorts.
- Adjust the CLI arguments (e.g., `--output-dir`) to keep separate audit trails for multiple studies.
- Embed the generated figure and summary table directly into assay validation documentation or slide decks to communicate findings.

Feel free to fork this repository, adapt the code to your own pipelines, and reach out if you want to collaborate on additional translational bioinformatics projects!
