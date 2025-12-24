"""Synthesis Agent for industrial synthesis route analysis."""
from typing import Dict, List
from .base_agent import BaseAgent


class SynthesisAgent(BaseAgent):
    """Agent specialized in chemical synthesis route planning."""

    def __init__(self):
        super().__init__(
            name="SynthesisAgent",
            role="Industrial synthesis route and raw materials specialist"
        )

    def get_synthesis_route(self, molecule_name: str) -> Dict[str, any]:
        """Get industrial synthesis route and raw materials.

        Args:
            molecule_name: Name of the molecule

        Returns:
            Dictionary with synthesis_route and raw_materials
        """
        system_instruction = (
            "You are an industrial chemistry process engineer. "
            "Provide accurate industrial-scale synthesis routes and raw materials. "
            "Always respond in valid JSON format only."
        )

        prompt = f"""Provide the industrial synthesis route for {molecule_name}.

Return ONLY valid JSON in this exact format:
{{
  "synthetic_route": [
    "Step 1: [Reaction description with reagents and conditions]",
    "Step 2: [Next step]",
    "..."
  ],
  "raw_materials": ["Material 1", "Material 2", "..."],
  "key_intermediates": ["Intermediate 1", "..."]
}}

Focus on the most common industrial process. Be specific and practical."""

        response = self._call_api(prompt, system_instruction)
        data = self._extract_json(response)

        if data:
            return data
        else:
            return {
                "synthetic_route": ["Synthesis route not available"],
                "raw_materials": [],
                "key_intermediates": [],
                "error": "Failed to parse synthesis data"
            }

    def estimate_complexity(self, synthesis_route: List[str]) -> str:
        """Estimate synthesis complexity.

        Args:
            synthesis_route: List of synthesis steps

        Returns:
            Complexity assessment string
        """
        num_steps = len(synthesis_route)

        if num_steps <= 2:
            return "Simple (1-2 steps)"
        elif num_steps <= 4:
            return "Moderate (3-4 steps)"
        elif num_steps <= 6:
            return "Complex (5-6 steps)"
        else:
            return f"Very Complex ({num_steps} steps)"

    def process(self, molecule_name: str) -> Dict:
        """Process complete synthesis analysis.

        Args:
            molecule_name: Name of the molecule

        Returns:
            Dictionary with synthesis analysis
        """
        synthesis_data = self.get_synthesis_route(molecule_name)
        complexity = self.estimate_complexity(
            synthesis_data.get("synthetic_route", [])
        )

        return {
            **synthesis_data,
            "complexity": complexity
        }
