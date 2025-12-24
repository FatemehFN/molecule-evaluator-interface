"""Chemistry Agent for molecular structure analysis."""
from typing import Dict, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from .base_agent import BaseAgent


class ChemistryAgent(BaseAgent):
    """Agent specialized in chemical structure analysis using RDKit and Gemini."""

    def __init__(self):
        super().__init__(
            name="ChemistryAgent",
            role="Molecular structure and property analysis specialist"
        )

    def name_to_smiles(self, molecule_name: str) -> Optional[str]:
        """Convert molecule name to SMILES string.

        Args:
            molecule_name: Name of the molecule

        Returns:
            SMILES string or None if not found
        """
        system_instruction = (
            "You are a chemistry database. Return ONLY the canonical SMILES string "
            "for the molecule requested. No explanations, no extra text."
        )

        prompt = f"SMILES for {molecule_name}:"
        smiles = self._call_api(prompt, system_instruction)

        # Validate SMILES
        if smiles and not smiles.startswith("Error"):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return smiles

        return None

    def get_molecular_properties(self, smiles: str) -> Tuple[Optional[Dict], Optional[str]]:
        """Get molecular properties from SMILES using RDKit.

        Args:
            smiles: SMILES string

        Returns:
            Tuple of (properties dict, error message)
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, "Invalid SMILES string"

        try:
            properties = {
                "Molecular Weight": round(Descriptors.MolWt(mol), 2),
                "Molecular Formula": rdMolDescriptors.CalcMolFormula(mol),
                "LogP": round(Descriptors.MolLogP(mol), 2),
                "H-Bond Donors": rdMolDescriptors.CalcNumHBD(mol),
                "H-Bond Acceptors": rdMolDescriptors.CalcNumHBA(mol),
                "Rotatable Bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                "TPSA": round(Descriptors.TPSA(mol), 2),
            }
            return properties, None
        except Exception as e:
            return None, str(e)

    def get_applications(self, molecule_name: str) -> str:
        """Get industrial applications for a molecule.

        Args:
            molecule_name: Name of the molecule

        Returns:
            Formatted applications text
        """
        system_instruction = (
            "You are a pharmaceutical and chemical industry expert. "
            "Provide concise, accurate applications for chemicals in "
            "medicine, supplements, cosmetics, and other industries."
        )

        prompt = f"""List the main industrial applications of {molecule_name}:

Categories to cover:
- Pharmaceuticals & Medicine
- Dietary Supplements
- Cosmetics & Skincare
- Industrial Uses
- Agriculture (if applicable)

Format as concise bullet points. Be specific and accurate."""

        return self._call_api(prompt, system_instruction)

    def process(self, molecule_name: str) -> Dict:
        """Process complete chemistry analysis for a molecule.

        Args:
            molecule_name: Name of the molecule

        Returns:
            Dictionary with chemistry analysis results
        """
        smiles = self.name_to_smiles(molecule_name)

        if not smiles:
            return {
                "name": molecule_name,
                "smiles": None,
                "properties": None,
                "applications": "Could not retrieve information",
                "error": "SMILES not found"
            }

        properties, error = self.get_molecular_properties(smiles)
        applications = self.get_applications(molecule_name)

        return {
            "name": molecule_name,
            "smiles": smiles,
            "properties": properties,
            "applications": applications,
            "error": error
        }
