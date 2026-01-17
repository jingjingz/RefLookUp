# 📚 RefLookUp

> **RefLookUp automates the tedious process of verifying bibliographies by validating raw references against PubMed and instantly filling in missing PMIDs, DOIs, and PMCIDs.**

Stop manually searching for identifiers. Turn messy or AI-generated citation lists into accurate, verified data with a single click.

## ✨ Features

- **Automated Verification**: Queries PubMed to find official PMIDs, PMCIDs, and DOIs.
- **Metadata Enrichment**: Fills in missing details like publication year, volume, issue, and journal names.
- **Format Standardization**: Convert references to **NLM** or **APA** style automatically.
- **Article Categorization**: Distinguishes between **Articles**, **Abstracts**, and **Preprints** (arXiv/biorxiv/medrxiv).
- **Format Adaptive**: Automatically detects and processes text, CSV, and BibTeX files.
- **Export Options**:
  - **CSV**: Spreadsheet-ready output with enriched metadata and ArXiv IDs.
  - **Text**: Formatted bibliographies (NLM/APA).
  - **RIS**: Direct import for EndNote, Zotero, and Mendeley.
- **ArXiv Support**: Automatically detects and formats arXiv preprints even if they aren't in PubMed.

## 🚀 Quick Start (Streamlit App)

The easiest way to use RefLookUp is via the web interface: **[https://reflookup.streamlit.app/](https://reflookup.streamlit.app/)**

### Running Locally

1. **Clone the repository:**

   ```bash
   git clone https://github.com/jingjingz/RefLookUp.git
   cd RefLookUp
   ```

2. **Install dependencies:**

   ```bash
   pip install -r requirements.txt
   ```

3. **Run the app:**

   ```bash
   ./run_app.sh
   # OR
   python3 -m streamlit run app.py
   ```

4. **Open your browser** to the URL shown (usually `http://localhost:8501`).

## 📖 Tutorial: Managing Citations for Biosketches

Use RefLookUp to verify your publications and easily import them into **MyNCBI** for your NIH Biosketch.

### 1. Export References

Export your selected references from Google Scholar (or another manager) as a `.bib` or `.csv` file.

![Google Scholar Export](Tutorial/googlescholar.png)

### 2. Verify in RefLookUp

Go to the **[Web App](https://reflookup.streamlit.app/)**, drag and drop your exported file into the upload field, and click **"Search References"**.
Once finished, look for the **RIS** download button. This format is specifically designed for citation managers.

![RefLookUp Interface](Tutorial/app.png)

### 3. Import to MyNCBI

Go to your **MyNCBI Bibliography**. Click **"Manage My Bibliography"** -> **"Add Citations"** -> **"From a file"**. Upload the `.ris` file you just downloaded.

![MyNCBI Import](Tutorial/myNCBI.png)

## 💻 CLI Usage (Batch Processing)

You can process multiple files at once directly from the terminal.

1. Place your input files (supported: `.txt`, `.csv`, `.bib`) inside the `input_examples/` folder.
2. Run the script:

   ```bash
   python3 ref_utils.py
   ```

3. The script will automatically process **all valid files** found in `input_examples/`.
4. Check the `output/` folder for your results (e.g., `citations_results.csv`, `mypapers_results.bib`).

## 🛠️ Project Structure

- `app.py`: The Streamlit web application.
- `ref_utils.py`: Core logic for PubMed searching, parsing, and formatting.
- `requirements.txt`: Python dependencies.
- `run_app.sh`: Helper script to launch the app reliably.

## 📄 License

MIT License. Free to use and modify.
