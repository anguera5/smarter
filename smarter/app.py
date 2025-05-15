import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit.Chem import Draw
import pandas as pd
from datetime import datetime

from utils import validate, generate_output, export_bytes

# Initialize session state
if "smarts" not in st.session_state:
    st.session_state.smarts = ""
if "compound_type" not in st.session_state:
    st.session_state.compound_type = "SMILES"

# Set configuration of the page
st.set_page_config(
    page_title="Smarter",
    page_icon="🪄",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'mailto:albert.anguera.sempere@gmail.com',
        'Report a bug': "https://github.com/anguera5/smarter/issues/new",
        'About': "**Smarter** is an interactive Streamlit web application designed to **simplify** the process of working with "
        "molecular substructure patterns. \n\nSmarter enables users to **input or draw SMARTS** patterns, "
        "validate them, and **apply them to lists of molecules** provided as SMILES or CAS numbers. The app **highlights substructure matches**, "
        "supports batch processing, and allows users to **export results to Excel**. "
        "\n\nSmarter is **ideal for chemists, researchers, and educators** who need a user-friendly tool for molecular pattern recognition and visualization."
    }
)

st.write("# Welcome to the Smarter App!")
st.write("## Generate a molecular pattern")
st.session_state.smarts = st.text_input("**Please enter your SMARTS pattern**", value=st.session_state.smarts, placeholder="SMARTS", help="Both regular [SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html) "
"and [ChemAxon's SMARTS](https://docs.chemaxon.com/display/docs/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.md) should work")

with st.expander("Or draw it", icon="✏️"):
    st.session_state.smarts = st_ketcher(st.session_state.smarts)
col1, col2 = st.columns([1,6], gap="small") 

with col1:
    smarts_button = st.button("Generate!", type="primary")

with col2:
    as_smiles = st.checkbox("Render as SMILES")

if smarts_button:
    pattern = validate(st.session_state.smarts, as_smiles=as_smiles)
    if pattern:
        st.image(Draw.MolToImage(pattern))
    else:
        st.error("Invalid SMARTS pattern. Please try again.")


st.write("## Apply the pattern to a set of molecules")

st.session_state.compound_type = st.radio(label="**Choose your input type**", options=["SMILES", "CAS"], horizontal=True,
                                          help="Use preferrably SMILES, since CAS needs to be transformed, which results in long computation times.")

input_placeholder = {"SMILES": "c1ccccc1\nCC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O\n...", "CAS": "7732-18-5\n80-05-7\n..."}
input_str = st.text_area(f"""**Please enter a {st.session_state.compound_type} (one {st.session_state.compound_type} per line) for 
                         which you would like to apply the pattern:**""", 
                         placeholder=input_placeholder[st.session_state.compound_type])
smiles_button = st.button("Apply!", type="primary")
if smiles_button:
    pattern = validate(st.session_state.smarts, as_smiles=False)
    if not pattern:
        st.error("Please enter a valid SMARTS pattern first.")
    output, invalid_smiles = generate_output(input_str, pattern, st.session_state.compound_type)
    frame = pd.DataFrame.from_dict(output, orient="index", columns=["Molecule"])
    frame.index.name = "SMILES"
    frame.reset_index(inplace=True)
    config = {
        "Molecule": st.column_config.ImageColumn(),
    }
    table = st.dataframe(frame, column_config=config, row_height=150)
    filename = f"smarter_pattern_result_{datetime.now().strftime("%Y_%m_%d_%H_%M_%S")}.xlsx"
    download_file = export_bytes(frame)
    
    if invalid_smiles:
        st.error(f"The following SMILES returned an error: {invalid_smiles}.")
    if table:
        st.download_button("Download Excel", download_file, filename, icon=":material/download:", mime="application/vnd.ms-excel")