from rdkit import Chem
from rdkit.Chem import Descriptors

def get_rdkit_info(molecule_name_or_smiles):
    mol = Chem.MolFromSmiles(molecule_name_or_smiles)
    if mol is None:
        return None, "Invalid molecule name or SMILES"

    return {
        "Name": molecule_name_or_smiles,
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "Molecular Formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "SMILES": molecule_name_or_smiles
    }, None
