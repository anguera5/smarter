from rdkit import Chem
from io import BytesIO
import base64
import pandas as pd
from rdkit.Chem import Draw, rdGeneralizedSubstruct

def validate(smarts, as_smiles=True):
    """
    Validate a given SMARTS or SMILES pattern.

    Parameters:
        smarts (str): The SMARTS or SMILES pattern to validate.
        as_smiles (bool, optional): If True, interpret the input as a SMILES pattern.
                                    If False, interpret it as a SMARTS pattern.
                                    Defaults to True.

    Returns:
        Union[rdkit.Chem.Mol, bool]: Returns an RDKit molecule object if the pattern is valid.
                                     Returns False if the pattern is invalid or an error occurs.
    """
    if not smarts:
        return False
    try:
        # Attempt to create a molecule from the SMARTS pattern
        mol = Chem.MolFromSmiles(smarts) if as_smiles else Chem.MolFromSmarts(smarts)
        if mol is None:
            return False
        return mol
    except Exception as e:
        return False

def generate_output(smiles_str, pattern):
    """
    Processes a list of SMILES strings, validates them, checks for substructure matches,
    and generates output with match results or visual representations.

    Args:
        smiles_str (str): A string containing multiple SMILES strings separated by newlines.
        pattern (str): A SMARTS pattern used for substructure matching.

    Returns:
        tuple: A tuple containing:
            - output (dict): A dictionary where keys are valid SMILES strings and values are either:
                - "No match" if the substructure pattern is not found.
                - A base64-encoded PNG image string of the molecule with the matched substructure highlighted.
            - invalid_smiles (list): A list of SMILES strings that are invalid.
    """
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
                bio = BytesIO()
                img.save(bio, format="PNG")
                img_bytes = bio.getvalue()
                base64_str = base64.b64encode(img_bytes).decode("utf-8")
                output[smiles] = f"data:image/png;base64,{base64_str}"
        else:
            invalid_smiles.append(smiles)
    return output, invalid_smiles
    
def export_bytes(df):
    """
    Exports a DataFrame to an Excel file in memory with embedded images.
    This function takes a pandas DataFrame containing SMILES strings and base64-encoded
    images, and writes them to an Excel file. The resulting Excel file is stored in a 
    BytesIO object, which can be used for further processing or downloading.
    Args:
        df (pandas.DataFrame): A DataFrame with the following columns:
            - 'SMILES': A string representing the SMILES notation of a molecule.
            - 'Molecule': A base64-encoded string representing an image of the molecule.
    Returns:
        BytesIO: A BytesIO object containing the Excel file data.
    Notes:
        - The Excel file is created using the xlsxwriter engine.
        - The 'Molecule' column is expected to contain base64-encoded images in the format
          "data:image/png;base64,<base64_data>".
        - If the 'Molecule' column does not contain valid base64 data, the cell will display "No match!".
        - Images are scaled to 50% of their original size when embedded in the Excel file.
    Raises:
        Exception: If an error occurs while decoding the base64 image or embedding it in the Excel file.
    """
    output = BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        workbook = writer.book
        worksheet = workbook.add_worksheet('Sheet1')
        writer.sheets['Sheet1'] = worksheet

        # Write header
        worksheet.write(0, 0, "SMILES")
        worksheet.write(0, 1, "Molecule")

        for idx, row in df.iterrows():
            worksheet.write(idx + 1, 0, row['SMILES'])
            worksheet.set_row(idx + 1, 100)

            # Decode base64 image
            try:
                header, b64data = row['Molecule'].split(',', 1)
                image_bytes = base64.b64decode(b64data)
                image_stream = BytesIO(image_bytes)

                # Save temporarily to insert (xlsxwriter needs file-like)
                worksheet.embed_image(idx + 1, 1, f'image_{idx}.png', {'image_data': image_stream, 'x_scale': 0.5, 'y_scale': 0.5})
            except Exception as e:
                worksheet.write(idx + 1, 1, "No match!")

        writer.close()
        return output