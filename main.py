import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def tanimoto_calc(smi1, smi2):
    """
    Calculates the Tanimoto similarity between two SMILES strings.

    Args:
        smi1 (str): The first SMILES string.
        smi2 (str): The second SMILES string.

    Returns:
        float: The Tanimoto similarity score between the two molecules, rounded to 3 decimal places.
    """

    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)

    if not mol1 or not mol2:
        return 0  # Handle cases where molecule conversion fails

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)

    try:
        s = round(DataStructs.TanimotoSimilarity(fp1, fp2), 3)
    except ZeroDivisionError:
        # Handle cases where both fingerprints are empty (no common bits)
        s = 0

    return s





 
