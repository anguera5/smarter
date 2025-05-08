from rdkit import Chem

def validate(smarts, as_smiles=True):
    """Validate the SMARTS pattern."""
    if not smarts:
        return False
    try:
        print(repr(smarts))
        # Attempt to create a molecule from the SMARTS pattern
        mol = Chem.MolFromSmiles(smarts) if as_smiles else Chem.MolFromSmarts(smarts)
        if mol is None:
            return False
        return mol
    except Exception as e:
        return False