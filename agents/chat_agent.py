"""Chat Agent for conversational interactions."""
from typing import List, Dict
from .base_agent import BaseAgent


class ChatAgent(BaseAgent):
    """Agent specialized in conversational interactions about molecules."""

    def __init__(self):
        super().__init__(
            name="ChatAgent",
            role="Conversational assistant for molecule analysis queries"
        )

    def chat(self, messages: List[Dict[str, str]], context: str = None) -> str:
        """Process chat conversation with context.

        Args:
            messages: List of message dictionaries with 'role' and 'content'
            context: Optional context about molecules

        Returns:
            Assistant response
        """
        # Build system instruction
        system_instruction = (
            "You are an expert chemistry and pharmaceutical industry assistant. "
            "You help analyze molecules, synthesis routes, market trends, and applications. "
            "Provide accurate, professional, and helpful responses."
        )

        if context:
            system_instruction += f"\n\nCONTEXT:\n{context}"

        # Build conversation prompt
        conversation = []
        for msg in messages:
            role = msg.get("role", "user")
            content = msg.get("content", "")

            if role == "user":
                conversation.append(f"User: {content}")
            elif role == "assistant":
                conversation.append(f"Assistant: {content}")

        conversation_text = "\n\n".join(conversation)
        conversation_text += "\n\nAssistant:"

        response = self._call_api(conversation_text, system_instruction)
        return response

    def format_molecule_context(self, molecule_data: List[Dict]) -> str:
        """Format molecule data for context.

        Args:
            molecule_data: List of molecule analysis results

        Returns:
            Formatted context string
        """
        if not molecule_data:
            return ""

        context_parts = []
        for mol in molecule_data:
            props = mol.get("properties", {})
            if not props:
                continue

            context = f"""
Molecule: {mol.get('name', 'Unknown')}
SMILES: {mol.get('smiles', 'N/A')}
Formula: {props.get('Molecular Formula', 'N/A')}
Molecular Weight: {props.get('Molecular Weight', 'N/A')} g/mol
LogP: {props.get('LogP', 'N/A')}

Synthesis Route:
{self._format_synthesis_route(mol.get('synthesis', {}))}

Raw Materials: {', '.join(mol.get('synthesis', {}).get('raw_materials', []))}

Price Info: {mol.get('market', {}).get('price_info', 'N/A')}
Total Cost: {mol.get('market', {}).get('total_price_toman', 0):,} Toman

Applications:
{mol.get('applications', 'N/A')}

Feasibility: {mol.get('feasibility', {}).get('feasibility_assessment', 'N/A')}
"""
            context_parts.append(context.strip())

        return "\n\n---\n\n".join(context_parts)

    def _format_synthesis_route(self, synthesis: Dict) -> str:
        """Format synthesis route for display.

        Args:
            synthesis: Synthesis data dictionary

        Returns:
            Formatted synthesis route
        """
        route = synthesis.get("synthetic_route", [])
        if not route:
            return "N/A"

        return "\n".join(route)

    def process(self, messages: List[Dict[str, str]], molecule_data: List[Dict] = None) -> str:
        """Process chat with molecule context.

        Args:
            messages: Conversation messages
            molecule_data: Optional molecule analysis data for context

        Returns:
            Chat response
        """
        context = None
        if molecule_data:
            context = self.format_molecule_context(molecule_data)

        return self.chat(messages, context)
