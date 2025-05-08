import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, rdGeneralizedSubstruct
import pandas as pd
import io
import base64

from utils import validate

st.write("# Welcome to the Smarter App!")
with st.expander("Please enter your SMARTS pattern:", icon="ℹ️"):
    st.info("**NOTE:**  Both regular SMARTS and ChemAxon's SMARTS should work")
smarts = st.text_input("text-input", placeholder="SMARTS", label_visibility="collapsed")
col1, col2 = st.columns([1,6], gap="small") 

with col1:
    smarts_button = st.button("Generate!", type="primary")

with col2:
    as_smiles = st.checkbox("Render as SMILES")

if smarts_button:
    pattern = validate(smarts, as_smiles=as_smiles)
    if pattern:
        st.image(Draw.MolToImage(pattern))
    else:
        st.error("Invalid SMARTS pattern. Please try again.")

smiles_str = st.text_area("Please enter a SMILES (one SMILES per line) for which you would like to apply the pattern:", 
                          placeholder="""c1ccccc1\nCC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O\n...""")
smiles_button = st.button("Apply!", type="primary")
if smiles_button:
    pattern = validate(smarts, as_smiles=False)
    if not pattern:
        st.error("Please enter a valid SMARTS pattern first.")

    invalid_smiles = []
    output = {}
    for smiles in smiles_str.split("\n"):
        mol = validate(smiles)
        if mol:
            xqm = rdGeneralizedSubstruct.CreateExtendedQueryMol(pattern)
            if not rdGeneralizedSubstruct.MolHasSubstructMatch(mol, xqm):
                output[smiles] = "No match"
            else:
                matches = rdGeneralizedSubstruct.MolGetSubstructMatches(mol, xqm)
                img = Draw.MolToImage(mol, highlightAtoms=matches[0])
                bio = io.BytesIO()
                img.save(bio, format="PNG")
                img_bytes = bio.getvalue()
                base64_str = base64.b64encode(img_bytes).decode("utf-8")
                output[smiles] = f"data:image/png;base64,{base64_str}"
        else:
            invalid_smiles.append(smiles)
    frame = pd.DataFrame.from_dict(output, orient="index", columns=["Molecule"])
    frame.index.name = "SMILES"
    frame.reset_index(inplace=True)
    config = {
        "Molecule": st.column_config.ImageColumn(),
    }
    st.dataframe(frame, column_config=config, row_height=150)
    if invalid_smiles:
        st.error(f"The following SMILES returned an error: {invalid_smiles}.")