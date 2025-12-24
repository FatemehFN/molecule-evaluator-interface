"""Fast Unified Agent with fallback mock data for marketing team testing."""
from typing import Dict, List
from .base_agent import BaseAgent
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import json


class UnifiedAgentFast(BaseAgent):
    """Ultra-fast agent with smart fallbacks for demo/testing."""

    def __init__(self):
        super().__init__(
            name="UnifiedAgentFast",
            role="Fast molecular analysis specialist"
        )
        self.mock_database = self._init_mock_database()

    def _init_mock_database(self) -> Dict:
        """Pre-loaded mock data for common molecules (for demo/fallback)."""
        return {
            "aspirin": {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "synthesis": {
                    "synthetic_route": [
                        "Step 1: Reaction of salicylic acid with acetic anhydride in presence of pyridine catalyst",
                        "Step 2: Crystallization and purification from acetic acid solution",
                    ],
                    "raw_materials": ["Salicylic acid", "Acetic anhydride", "Pyridine"],
                    "key_intermediates": ["Salicylic acid"],
                    "complexity": "Simple (2 steps)"
                },
                "applications": {
                    "pharmaceuticals": "Pain relief, anti-inflammatory, antiplatelet therapy",
                    "supplements": "General wellness, cardiovascular support",
                    "cosmetics": "Exfoliating agents",
                    "industrial": "Chemical synthesis intermediate"
                },
                "market": {
                    "raw_material_prices_toman": [
                        "Salicylic acid: 2,000,000-3,000,000 Toman/kg",
                        "Acetic anhydride: 5,000,000-7,000,000 Toman/kg"
                    ],
                    "estimated_total_cost_toman": 4500000
                }
            },
            "vitamin c": {
                "smiles": "O=C(O)[C@H](O)[C@H](O)[C@@H](O)CO",
                "synthesis": {
                    "synthetic_route": [
                        "Step 1: Glucose fermentation to gluconic acid",
                        "Step 2: Lactonization to glucono-delta-lactone",
                        "Step 3: Sorbitol reduction",
                        "Step 4: Oxidation to L-ascorbic acid",
                    ],
                    "raw_materials": ["Glucose", "Sorbitol", "Acetone", "Catalysts"],
                    "key_intermediates": ["Glucono-delta-lactone", "2-Keto-L-gulonic acid"],
                    "complexity": "Moderate (4 steps)"
                },
                "applications": {
                    "pharmaceuticals": "Immune support, antioxidant, collagen synthesis",
                    "supplements": "Dietary supplements, multivitamins",
                    "cosmetics": "Anti-aging, brightening agents",
                    "industrial": "Food preservation, pharmaceutical intermediate"
                },
                "market": {
                    "raw_material_prices_toman": [
                        "Glucose: 3,000,000-4,500,000 Toman/ton",
                        "Sorbitol: 8,000,000-10,000,000 Toman/ton"
                    ],
                    "estimated_total_cost_toman": 6000000
                }
            },
            "retinol": {
                "smiles": "CC(C)C=CC(C)=CC=CC(C)=CC=CC=C(C)C=CC(C)=CC(C)(C)C",
                "synthesis": {
                    "synthetic_route": [
                        "Step 1: Wittig reaction of ionone with triphenylphosphine alkyl",
                        "Step 2: Isomerization to trans-retinol",
                        "Step 3: Purification and stabilization",
                    ],
                    "raw_materials": ["Beta-ionone", "Phosphorus reagents", "Solvents"],
                    "key_intermediates": ["Retinal", "Retinoic acid"],
                    "complexity": "Complex (3 steps)"
                },
                "applications": {
                    "pharmaceuticals": "Anti-aging, skin regeneration, acne treatment",
                    "supplements": "Vision support, immune function",
                    "cosmetics": "Wrinkle reduction, skin texture improvement",
                    "industrial": "Cosmetic ingredient, pharmaceutical active"
                },
                "market": {
                    "raw_material_prices_toman": [
                        "Beta-ionone: 25,000,000-35,000,000 Toman/kg",
                        "Phosphorus compounds: 15,000,000-20,000,000 Toman/kg"
                    ],
                    "estimated_total_cost_toman": 28000000
                }
            },
            "niacinamide": {
                "smiles": "c1cc(cnc1)C(=O)N",
                "synthesis": {
                    "synthetic_route": [
                        "Step 1: Oxidation of nicotine to nicotinic acid",
                        "Step 2: Amidation with ammonia",
                        "Step 3: Crystallization and purification",
                    ],
                    "raw_materials": ["Nicotinic acid", "Ammonia", "Catalysts"],
                    "key_intermediates": ["Nicotinic acid"],
                    "complexity": "Simple (3 steps)"
                },
                "applications": {
                    "pharmaceuticals": "B vitamin supplementation, skin barrier support",
                    "supplements": "Energy metabolism, general wellness",
                    "cosmetics": "Pore refinement, brightening, anti-inflammatory",
                    "industrial": "Food fortification, cosmetic ingredient"
                },
                "market": {
                    "raw_material_prices_toman": [
                        "Nicotinic acid: 5,000,000-8,000,000 Toman/kg",
                        "Ammonia: 2,000,000-3,000,000 Toman/1000L"
                    ],
                    "estimated_total_cost_toman": 5500000
                }
            }
        }

    def get_molecular_properties_local(self, smiles: str) -> Dict:
        """Calculate molecular properties locally with RDKit."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}

        try:
            return {
                "Molecular Weight": round(Descriptors.MolWt(mol), 2),
                "Molecular Formula": rdMolDescriptors.CalcMolFormula(mol),
                "LogP": round(Descriptors.MolLogP(mol), 2),
                "H-Bond Donors": rdMolDescriptors.CalcNumHBD(mol),
                "H-Bond Acceptors": rdMolDescriptors.CalcNumHBA(mol),
                "Rotatable Bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                "TPSA": round(Descriptors.TPSA(mol), 2),
            }
        except Exception as e:
            return {"error": str(e)}

    def process(self, molecule_name: str) -> Dict:
        """Fast analysis using REST API with instant fallback for known molecules."""
        import sys
        
        mol_key = molecule_name.lower().strip()
        
        # Check cache first for instant response on known molecules
        if mol_key in self.mock_database:
            print(f"⚡ Using cached data for {molecule_name}", file=sys.stderr, flush=True)
            data = self.mock_database[mol_key]
        else:
            # Use REST API for unknown molecules
            print(f"🔍 Fetching data for {molecule_name} from API...", file=sys.stderr, flush=True)
            data = self._fetch_from_api(molecule_name)

        if not data or not data.get("smiles"):
            return {
                "name": molecule_name,
                "status": "failed",
                "error": "Could not retrieve molecule data"
            }

        # Calculate properties locally
        smiles = data.get("smiles")
        properties = self.get_molecular_properties_local(smiles)

        # Format synthesis
        synthesis = data.get("synthesis", {})
        
        # Format applications
        applications = data.get("applications", {})
        apps_text = []
        if applications.get("pharmaceuticals"):
            apps_text.append(f"**Pharmaceuticals:** {applications['pharmaceuticals']}")
        if applications.get("supplements"):
            apps_text.append(f"**Supplements:** {applications['supplements']}")
        if applications.get("cosmetics"):
            apps_text.append(f"**Cosmetics:** {applications['cosmetics']}")
        if applications.get("industrial"):
            apps_text.append(f"**Industrial:** {applications['industrial']}")
        applications_formatted = "\n\n".join(apps_text) if apps_text else "No application data"

        # Format market
        market_data = data.get("market", {})
        price_list = market_data.get("raw_material_prices_toman", [])
        price_info = "\n".join(f"• {p}" for p in price_list) if price_list else "Price info unavailable"
        total_price = market_data.get("estimated_total_cost_toman", 0)

        market = {
            "price_info": price_info,
            "total_price_toman": total_price
        }

        # Quick feasibility
        feasibility = {
            "feasibility_assessment": "✅ Production feasible with standard equipment",
            "equipment_available": True
        }

        return {
            "name": molecule_name,
            "smiles": smiles,
            "properties": properties,
            "synthesis": synthesis,
            "market": market,
            "feasibility": feasibility,
            "applications": applications_formatted,
            "status": "success"
        }

    def _fetch_from_api(self, molecule_name: str) -> Dict:
        """Call Gemini API via REST for molecule data."""
        system_instruction = (
            "You are a chemistry expert. Provide accurate molecular data in valid JSON only. "
            "No markdown, no explanations, just JSON."
        )

        prompt = f"""For the molecule "{molecule_name}", provide this exact JSON structure:
{{
  "smiles": "canonical SMILES string",
  "synthesis": {{
    "synthetic_route": ["Step 1 with reagents", "Step 2"],
    "raw_materials": ["Material 1", "Material 2"],
    "key_intermediates": ["Intermediate 1"],
    "complexity": "Simple/Moderate/Complex/Very Complex"
  }},
  "applications": {{
    "pharmaceuticals": "Medical uses",
    "supplements": "Supplement uses",
    "cosmetics": "Cosmetic uses",
    "industrial": "Industrial uses"
  }},
  "market": {{
    "raw_material_prices_toman": ["Item: price range"],
    "estimated_total_cost_toman": 0
  }}
}}

Return ONLY valid JSON. No markdown. No extra text."""

        response = self._call_api(prompt, system_instruction, timeout=60)
        return self._extract_json(response)
