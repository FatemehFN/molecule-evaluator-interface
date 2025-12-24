"""Unified Agent - Combines multiple analyses in single API calls for speed."""
from typing import Dict, List
from .base_agent import BaseAgent
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


class UnifiedAgent(BaseAgent):
    """Fast agent that combines multiple analyses in fewer API calls."""

    def __init__(self):
        super().__init__(
            name="UnifiedAgent",
            role="Unified molecular analysis specialist"
        )

    def get_molecular_properties_local(self, smiles: str) -> Dict:
        """Calculate molecular properties locally with RDKit (no API call).

        Args:
            smiles: SMILES string

        Returns:
            Dictionary with properties or error
        """
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

    def analyze_complete(self, molecule_name: str) -> Dict:
        """Complete analysis in just 2 API calls instead of 5.

        Call 1: Get SMILES + Synthesis + Applications (combined)
        Call 2: Get Feasibility (needs synthesis results)
        Local: Calculate properties with RDKit (no API needed)

        Args:
            molecule_name: Name of the molecule

        Returns:
            Complete analysis dictionary
        """
        import sys
        print(f"🔍 Processing: {molecule_name}", file=sys.stderr, flush=True)
        
        # CALL 1: Combined SMILES, Synthesis, and Applications
        system_instruction = (
            "You are a chemistry expert. Provide accurate molecular data, "
            "industrial synthesis routes, and applications. Always return valid JSON."
        )

        prompt = f"""For the molecule **{molecule_name}**, provide a comprehensive analysis in JSON format:

{{
  "smiles": "canonical SMILES string",
  "synthesis": {{
    "synthetic_route": [
      "Step 1: [Detailed reaction with reagents and conditions]",
      "Step 2: [Next step]",
      "..."
    ],
    "raw_materials": ["Material 1", "Material 2", "..."],
    "key_intermediates": ["Intermediate 1", "..."],
    "complexity": "Simple/Moderate/Complex/Very Complex"
  }},
  "applications": {{
    "pharmaceuticals": "Uses in medicine/pharma",
    "supplements": "Uses in dietary supplements",
    "cosmetics": "Uses in skincare/cosmetics",
    "industrial": "Other industrial uses"
  }},
  "market": {{
    "raw_material_prices_toman": ["Material: price range", "..."],
    "estimated_total_cost_toman": 0
  }}
}}

Return ONLY the JSON, no markdown or extra text."""

        response = self._call_api(prompt, system_instruction)
        import sys
        print(f"📤 Raw API Response (first 500 chars): {response[:500]}", file=sys.stderr, flush=True)
        
        data = self._extract_json(response)

        if not data or not data.get("smiles"):
            return {
                "name": molecule_name,
                "status": "failed",
                "error": "Could not retrieve molecule data",
                "raw_response": response[:200]
            }

        smiles = data.get("smiles")

        # Calculate properties locally (no API call needed!)
        properties = self.get_molecular_properties_local(smiles)

        # Format synthesis data
        synthesis = data.get("synthesis", {})
        complexity = synthesis.get("complexity", "Unknown")
        if not complexity or complexity == "Unknown":
            num_steps = len(synthesis.get("synthetic_route", []))
            if num_steps <= 2:
                complexity = "Simple (1-2 steps)"
            elif num_steps <= 4:
                complexity = "Moderate (3-4 steps)"
            elif num_steps <= 6:
                complexity = "Complex (5-6 steps)"
            else:
                complexity = f"Very Complex ({num_steps} steps)"
        synthesis["complexity"] = complexity

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

        applications_formatted = "\n\n".join(apps_text) if apps_text else "No application data available"

        # Format market data
        market_data = data.get("market", {})
        price_list = market_data.get("raw_material_prices_toman", [])
        price_info = "\n".join(f"• {p}" for p in price_list) if price_list else "Price information not available"

        # Calculate total from price ranges
        total_price = self._calculate_price_from_text(price_info)

        market = {
            "price_info": price_info,
            "total_price_toman": total_price
        }

        # CALL 2: Feasibility (needs synthesis route)
        feasibility = self._assess_feasibility_fast(
            synthesis.get("synthetic_route", [])
        )

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

    def _assess_feasibility_fast(self, synthesis_steps: List[str]) -> Dict:
        """Quick feasibility assessment.

        Args:
            synthesis_steps: List of synthesis steps

        Returns:
            Feasibility dictionary
        """
        if not synthesis_steps or synthesis_steps == ["Synthesis route not available"]:
            return {
                "feasibility_assessment": "❌ Cannot assess - synthesis route unavailable",
                "equipment_available": False
            }

        # Load equipment list
        from pathlib import Path
        from config import Config

        equipment_path = Path(Config.EQUIPMENT_LIST_PATH)
        equipment_list = []

        if equipment_path.exists():
            try:
                with open(equipment_path, "r", encoding="utf-8") as f:
                    equipment_list = [line.strip("-• \n") for line in f if line.strip()]
            except Exception:
                pass

        if not equipment_list:
            return {
                "feasibility_assessment": "⚠️ Equipment list not available",
                "equipment_available": False
            }

        # Simple feasibility check
        equipment_text = "\n".join(f"- {eq}" for eq in equipment_list[:10])  # Just first 10
        route_text = "\n".join(f"{i+1}. {step}" for i, step in enumerate(synthesis_steps[:5]))  # Just first 5 steps

        system_instruction = (
            "You are a chemical production analyst. Assess feasibility concisely."
        )

        prompt = f"""Quick feasibility assessment:

SYNTHESIS ROUTE:
{route_text}

EQUIPMENT AVAILABLE:
{equipment_text}

Reply in ONE sentence with:
✅ FEASIBLE if compatible
⚠️ PARTIALLY FEASIBLE if some equipment missing
❌ NOT FEASIBLE if critical equipment missing

Be concise."""

        assessment = self._call_api(prompt, system_instruction)

        return {
            "feasibility_assessment": assessment.strip(),
            "equipment_available": True
        }

    def _calculate_price_from_text(self, price_text: str) -> int:
        """Calculate total average price from Persian price text.

        Args:
            price_text: Price information

        Returns:
            Total average price in Toman
        """
        import re

        # Convert Persian digits to English
        persian_to_english = str.maketrans("۰۱۲۳۴۵۶۷۸۹٬", "0123456789,")
        normalized = price_text.translate(persian_to_english)

        total = 0
        matches = re.findall(r"(\d[\d,]*)\s*[-–—]\s*(\d[\d,]*)", normalized)

        for low_str, high_str in matches:
            try:
                low = int(low_str.replace(",", ""))
                high = int(high_str.replace(",", ""))
                total += (low + high) / 2
            except ValueError:
                continue

        return round(total)

    def process(self, molecule_name: str) -> Dict:
        """Process molecule with optimized speed.

        Args:
            molecule_name: Name of the molecule

        Returns:
            Complete analysis
        """
        return self.analyze_complete(molecule_name)
