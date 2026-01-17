# 📚 RefLookUp

> **RefLookUp automates the tedious process of verifying bibliographies by validating raw references against PubMed and instantly filling in missing PMIDs, DOIs, and PMCIDs.**

Stop manually searching for identifiers. Turn messy or AI-generated citation lists into accurate, verified data with a single click.

## ✨ Features

- **Automated Verification**: Queries PubMed to find official PMIDs, PMCIDs, and DOIs.
- **Metadata Enrichment**: Fills in missing details like publication year, volume, issue, and journal names.
- **Format Standardization**: Convert references to **NLM** or **APA** style automatically.
- **Article Categorization**: Distinguishes between full **Articles** and **Abstracts/Posters** based on publication type.
- **Export Options**:
  - **CSV**: For spreadsheet analysis (includes ArXiv IDs).
  - **Text**: Formatted bibliography ready for word processors.
  - **RIS**: Import directly into EndNote, Zotero, or Mendeley.
- **ArXiv Support**: Automatically detects and formats arXiv preprints even if they aren't in PubMed.

## 🚀 Quick Start (Streamlit App)

The easiest way to use RefLookUp is via the web interface.

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

### ☁️ Deploying to Streamlit Cloud

1. Push this repository to your GitHub.
2. Go to [share.streamlit.io](https://share.streamlit.io/).
3. Create a **New App**.
4. Select this repository and set the main file to `app.py`.
5. **Important**: When using the hosted app, you MUST provide your email address in the sidebar so PubMed doesn't rate-limit the shared server IP.

## 💻 CLI Usage

You can also run the tool directly from the terminal for batch processing.

1. Create a file named `references.txt` with your list of citations.
2. Run the script:

   ```bash
   python3 refcheck.py
   ```

3. Check the output files:
   - `recent_references.csv`
   - `recent_references.txt`
   - `recent_references.ris`

## 🛠️ Project Structure

- `app.py`: The Streamlit web application.
- `refcheck.py`: Core logic for PubMed searching, parsing, and formatting.
- `requirements.txt`: Python dependencies.
- `run_app.sh`: Helper script to launch the app reliably.

## 📄 License

MIT License. Free to use and modify.
