import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit.Chem import Draw
import pandas as pd
from datetime import datetime

from utils import validate, generate_output, export_bytes

st.write("# Welcome to the Smarter App!")
if "smarts" not in st.session_state:
    st.session_state.smarts = ""
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

smiles_str = st.text_area("Please enter a SMILES (one SMILES per line) for which you would like to apply the pattern:", 
                          placeholder="""c1ccccc1\nCC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O\n...""")
smiles_button = st.button("Apply!", type="primary")
if smiles_button:
    pattern = validate(st.session_state.smarts, as_smiles=False)
    if not pattern:
        st.error("Please enter a valid SMARTS pattern first.")

    output, invalid_smiles = generate_output(smiles_str, pattern)
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