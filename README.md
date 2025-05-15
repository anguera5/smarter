# Smarter

**Smarter** is a Streamlit web application for rendering and applying SMARTS patterns to molecules using RDKit. It allows users to input or draw SMARTS patterns, validate them, and apply them to sets of molecules (SMILES or CAS numbers). The app highlights substructure matches and enables exporting results (including images) to Excel.

## Features

- **SMARTS Pattern Input:** Enter or draw SMARTS patterns interactively.
- **Validation:** Validate SMARTS or SMILES patterns using RDKit.
- **Batch Processing:** Apply patterns to lists of SMILES or CAS numbers.
- **Substructure Highlighting:** Visualize substructure matches in molecules.
- **Export:** Download results as an Excel file with embedded molecule images.
- **CAS to SMILES Conversion:** Automatically resolves CAS numbers to SMILES using Cirpy.

## Getting Started

### Prerequisites

- Python 3.12+
- [RDKit](https://www.rdkit.org/)
- [Streamlit](https://streamlit.io/)
- [Cirpy](https://github.com/mcs07/Cirpy)
- [Pandas](https://pandas.pydata.org/)
- [XlsxWriter](https://xlsxwriter.readthedocs.io/)
- [streamlit-ketcher](https://github.com/epam/streamlit-ketcher)

All dependencies are listed in [requirements.txt](requirements.txt) and [pyproject.toml](pyproject.toml).

### Installation

1. **Clone the repository:**
   ```sh
   git clone <your-repo-url>
   cd smarter
   ```

2. **Install system dependencies (if using devcontainer):**
   ```sh
   sudo apt update && sudo apt install -y libxrender1
   ```

3. **Install Python dependencies:**
   ```sh
   pip install -r requirements.txt
   ```

   Or, if using Poetry:
   ```sh
   poetry install
   ```

### Running the App

```sh
streamlit run smarter/app.py
```

The app will be available at [http://localhost:8501](http://localhost:8501).

## Usage

1. **Enter or draw a SMARTS pattern** in the input box or using the Ketcher widget.
2. **Validate and visualize** the pattern.
3. **Choose input type** (SMILES or CAS) and paste your list (one per line).
4. **Apply the pattern** to see which molecules match and visualize the results.
5. **Download the results** as an Excel file with images.

## Project Structure

- [app.py](smarter/app.py): Main Streamlit application.
- [utils.py](smarter/utils.py): Utility functions for validation, processing, and exporting.
- [requirements.txt](./requirements.txt): Python dependencies.
- [pyproject.toml](./pyproject.toml): Project metadata and dependencies.
- [devcontainer.json](.devcontainer/devcontainer.json): Devcontainer configuration for VS Code.
- [config.toml](.streamlit/config.toml): Streamlit theme configuration.

## Development

- The app is designed for use in VS Code with Dev Containers.
- On Codespaces or in a devcontainer, dependencies and the app will start automatically.

## License

MIT License

---

**Author:** Albert Anguera Sempere  
For questions, contact [albert.anguera.sempere@gmail.com](mailto:albert.anguera.sempere@gmail.com)
