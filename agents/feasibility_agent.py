"""Feasibility Agent for production feasibility analysis."""
from typing import List, Dict
from pathlib import Path
from .base_agent import BaseAgent
from config import Config


class FeasibilityAgent(BaseAgent):
    """Agent specialized in production feasibility assessment."""

    def __init__(self):
        super().__init__(
            name="FeasibilityAgent",
            role="Production feasibility and equipment compatibility specialist"
        )
        self.equipment_list = self._load_equipment_list()

    def _load_equipment_list(self) -> List[str]:
        """Load equipment list from file.

        Returns:
            List of equipment items
        """
        equipment_path = Path(Config.EQUIPMENT_LIST_PATH)

        if not equipment_path.exists():
            return []

        try:
            with open(equipment_path, "r", encoding="utf-8") as f:
                lines = [line.strip("-• \n") for line in f if line.strip()]
            return lines
        except Exception:
            return []

    def assess_feasibility(
        self,
        synthesis_steps: List[str],
        equipment_list: List[str] = None
    ) -> str:
        """Assess production feasibility based on synthesis and equipment.

        Args:
            synthesis_steps: List of synthesis route steps
            equipment_list: Optional custom equipment list

        Returns:
            Feasibility assessment text
        """
        equipment = equipment_list if equipment_list else self.equipment_list

        if not equipment:
            return "⚠️ Equipment list not available for assessment"

        if not synthesis_steps or synthesis_steps == ["Synthesis route not available"]:
            return "❌ Cannot assess feasibility - synthesis route not available"

        equipment_text = "\n".join(f"- {eq}" for eq in equipment)
        route_text = "\n".join(f"{i+1}. {step}" for i, step in enumerate(synthesis_steps))

        system_instruction = (
            "You are a chemical production feasibility analyst. "
            "Assess whether synthesis routes can be executed with given equipment. "
            "Be practical and specific."
        )

        prompt = f"""Assess production feasibility for this synthesis route:

SYNTHESIS ROUTE:
{route_text}

AVAILABLE EQUIPMENT:
{equipment_text}

Provide assessment in this format:

✅ FEASIBLE - [Brief explanation why it's feasible]

OR

⚠️ PARTIALLY FEASIBLE - [What's available, what's missing]

OR

❌ NOT FEASIBLE - [What critical equipment/capabilities are missing]

Be specific about equipment compatibility. Consider:
- Reaction vessels and reactors
- Temperature control
- Mixing and agitation
- Separation equipment
- Safety equipment"""

        return self._call_api(prompt, system_instruction)

    def get_recommendations(self, synthesis_steps: List[str]) -> str:
        """Get recommendations for production optimization.

        Args:
            synthesis_steps: List of synthesis steps

        Returns:
            Recommendations text
        """
        route_text = "\n".join(f"{i+1}. {step}" for i, step in enumerate(synthesis_steps))

        system_instruction = (
            "You are a chemical process optimization expert. "
            "Provide practical recommendations for industrial production."
        )

        prompt = f"""For this synthesis route, provide optimization recommendations:

{route_text}

Provide recommendations for:
1. Process efficiency improvements
2. Cost reduction opportunities
3. Safety considerations
4. Quality control points

Keep recommendations practical and industry-focused."""

        return self._call_api(prompt, system_instruction)

    def process(self, synthesis_steps: List[str], get_recommendations: bool = False) -> Dict:
        """Process complete feasibility analysis.

        Args:
            synthesis_steps: List of synthesis steps
            get_recommendations: Whether to include recommendations

        Returns:
            Dictionary with feasibility analysis
        """
        assessment = self.assess_feasibility(synthesis_steps)

        result = {
            "feasibility_assessment": assessment,
            "equipment_available": len(self.equipment_list) > 0
        }

        if get_recommendations:
            result["recommendations"] = self.get_recommendations(synthesis_steps)

        return result
