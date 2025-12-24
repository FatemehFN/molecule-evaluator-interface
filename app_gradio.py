"""MoleculeAI - Gradio Interface for Internal Use.

A fast, beautiful, and intuitive chatbot for analyzing molecules in chemical
ingredient companies. Built with Gradio for superior UX and performance.
"""

import gradio as gr
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
from typing import List, Dict, Tuple, Optional

from config import Config
from agents import (
    ChemistryAgent,
    SynthesisAgent,
    MarketAgent,
    FeasibilityAgent,
    ChatAgent
)


# Initialize configuration
try:
    Config.validate()
except ValueError as e:
    print(f"⚠️ Configuration Error: {e}")
    print("Please set GEMINI_API_KEY in your .env file")
    exit(1)


# Initialize agents (singleton pattern)
class AgentManager:
    """Manages all AI agents as singletons."""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.chemistry_agent = ChemistryAgent()
            cls._instance.synthesis_agent = SynthesisAgent()
            cls._instance.market_agent = MarketAgent()
            cls._instance.feasibility_agent = FeasibilityAgent()
            cls._instance.chat_agent = ChatAgent()
            cls._instance.molecule_data = []
        return cls._instance


agents = AgentManager()


def process_molecule(molecule_name: str) -> Dict:
    """Process a single molecule through all agents.

    Args:
        molecule_name: Name of the molecule

    Returns:
        Dictionary with complete analysis
    """
    result = {"name": molecule_name}

    # Chemistry analysis
    chemistry_data = agents.chemistry_agent.process(molecule_name)
    result.update(chemistry_data)

    if not chemistry_data.get("smiles"):
        result["status"] = "failed"
        return result

    # Synthesis analysis
    synthesis_data = agents.synthesis_agent.process(molecule_name)
    result["synthesis"] = synthesis_data

    # Market analysis
    market_data = agents.market_agent.process(
        raw_materials=synthesis_data.get("raw_materials", []),
        molecule_name=molecule_name
    )
    result["market"] = market_data

    # Feasibility analysis
    feasibility_data = agents.feasibility_agent.process(
        synthesis_steps=synthesis_data.get("synthetic_route", [])
    )
    result["feasibility"] = feasibility_data

    result["status"] = "success"
    return result


def analyze_molecules(file) -> Tuple[str, str, pd.DataFrame]:
    """Analyze molecules from uploaded file.

    Args:
        file: Uploaded text file

    Returns:
        Tuple of (status message, detailed results, dataframe)
    """
    if file is None:
        return "⚠️ Please upload a file first", "", pd.DataFrame()

    # Read molecule names
    content = file.decode('utf-8') if isinstance(file, bytes) else file
    if hasattr(file, 'name'):
        with open(file.name, 'r') as f:
            content = f.read()

    molecule_names = [line.strip() for line in content.split('\n') if line.strip()]

    if not molecule_names:
        return "⚠️ No molecules found in file", "", pd.DataFrame()

    # Process molecules
    results = []
    status_msg = f"🔬 Processing {len(molecule_names)} molecules...\n\n"

    for idx, name in enumerate(molecule_names, 1):
        status_msg += f"[{idx}/{len(molecule_names)}] Analyzing {name}...\n"
        result = process_molecule(name)
        results.append(result)

    agents.molecule_data = results

    # Create summary
    successful = sum(1 for r in results if r.get("status") == "success")
    status_msg += f"\n✅ Complete! {successful}/{len(molecule_names)} molecules analyzed successfully."

    # Create detailed results text
    detailed = format_detailed_results(results)

    # Create dataframe for table view
    df = create_summary_dataframe(results)

    return status_msg, detailed, df


def format_detailed_results(results: List[Dict]) -> str:
    """Format detailed analysis results.

    Args:
        results: List of molecule analysis results

    Returns:
        Formatted markdown text
    """
    output = "# 📊 Detailed Analysis Results\n\n"

    for mol in results:
        output += f"## 🧪 {mol.get('name', 'Unknown')}\n\n"

        if mol.get("status") == "failed":
            output += "❌ Analysis failed - molecule not found\n\n"
            output += "---\n\n"
            continue

        # Basic Properties
        props = mol.get("properties", {})
        output += "### Chemical Properties\n"
        output += f"- **SMILES:** `{mol.get('smiles', 'N/A')}`\n"
        output += f"- **Formula:** {props.get('Molecular Formula', 'N/A')}\n"
        output += f"- **Molecular Weight:** {props.get('Molecular Weight', 'N/A')} g/mol\n"
        output += f"- **LogP:** {props.get('LogP', 'N/A')}\n"
        output += f"- **H-Bond Donors:** {props.get('H-Bond Donors', 'N/A')}\n"
        output += f"- **H-Bond Acceptors:** {props.get('H-Bond Acceptors', 'N/A')}\n\n"

        # Synthesis
        synthesis = mol.get("synthesis", {})
        output += "### Synthesis Information\n"
        output += f"- **Complexity:** {synthesis.get('complexity', 'N/A')}\n"
        output += "- **Synthesis Route:**\n"
        for step in synthesis.get("synthetic_route", []):
            output += f"  - {step}\n"
        output += f"- **Raw Materials:** {', '.join(synthesis.get('raw_materials', []))}\n\n"

        # Market
        market = mol.get("market", {})
        total_price = market.get("total_price_toman", 0)
        if total_price > 0:
            output += f"### Market Information\n"
            output += f"- **Estimated Cost:** {total_price:,} Toman\n\n"

        # Feasibility
        feasibility = mol.get("feasibility", {})
        assessment = feasibility.get("feasibility_assessment", "")
        if assessment:
            output += f"### Production Feasibility\n{assessment}\n\n"

        # Applications
        applications = mol.get("applications", "")
        if applications and not applications.startswith("Could not"):
            output += f"### Applications\n{applications}\n\n"

        output += "---\n\n"

    return output


def create_summary_dataframe(results: List[Dict]) -> pd.DataFrame:
    """Create summary dataframe for table view.

    Args:
        results: List of molecule analysis results

    Returns:
        Pandas DataFrame
    """
    data = []

    for mol in results:
        if mol.get("status") == "failed":
            data.append({
                "Molecule": mol.get("name", "Unknown"),
                "Formula": "N/A",
                "MW (g/mol)": "N/A",
                "Complexity": "N/A",
                "Cost (Toman)": "N/A",
                "Status": "❌ Failed"
            })
        else:
            props = mol.get("properties", {})
            synthesis = mol.get("synthesis", {})
            market = mol.get("market", {})

            data.append({
                "Molecule": mol.get("name", "Unknown"),
                "Formula": props.get("Molecular Formula", "N/A"),
                "MW (g/mol)": props.get("Molecular Weight", "N/A"),
                "Complexity": synthesis.get("complexity", "N/A"),
                "Cost (Toman)": f"{market.get('total_price_toman', 0):,}",
                "Status": "✅ Success"
            })

    return pd.DataFrame(data)


def chat_with_ai(message: str, history: List) -> Tuple[str, List]:
    """Chat with AI assistant about molecules.

    Args:
        message: User message
        history: Chat history

    Returns:
        Tuple of (empty string to clear input, updated history)
    """
    if not message.strip():
        return "", history

    if not agents.molecule_data:
        response = "Please analyze some molecules first before asking questions."
        history.append((message, response))
        return "", history

    # Build message list for chat agent
    messages = []
    for user_msg, bot_msg in history:
        messages.append({"role": "user", "content": user_msg})
        messages.append({"role": "assistant", "content": bot_msg})
    messages.append({"role": "user", "content": message})

    # Get AI response
    response = agents.chat_agent.process(messages, agents.molecule_data)

    # Update history
    history.append((message, response))

    return "", history


def create_market_chart(molecule_name: str) -> Optional[plt.Figure]:
    """Create market trend chart for a molecule.

    Args:
        molecule_name: Name of the molecule

    Returns:
        Matplotlib figure or None
    """
    if not molecule_name:
        return None

    # Get market trend data
    trend_data = agents.market_agent.get_market_trend(molecule_name)

    if not trend_data or not trend_data.get("global_market"):
        return None

    # Create chart
    fig, ax = plt.subplots(figsize=(10, 6))

    years = sorted(trend_data["global_market"].keys())
    global_vals = [trend_data["global_market"][y] for y in years]
    me_vals = [trend_data["middle_east_market"][y] for y in years]

    ax.plot(years, global_vals, marker='o', linewidth=3,
            label="Global Market", color="#667eea", markersize=8)
    ax.plot(years, me_vals, marker='s', linewidth=3,
            label="Middle East Market", color="#764ba2",
            linestyle="--", markersize=8)

    ax.set_title(f"Market Value Projection: {molecule_name}",
                fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel("Year", fontsize=13)
    ax.set_ylabel("Market Value (Million USD)", fontsize=13)
    ax.legend(loc='best', fontsize=12, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()
    return fig


def create_cost_comparison() -> Optional[plt.Figure]:
    """Create cost comparison chart.

    Returns:
        Matplotlib figure or None
    """
    if not agents.molecule_data:
        return None

    # Get molecules with costs
    molecules_with_costs = [
        m for m in agents.molecule_data
        if m.get("market", {}).get("total_price_toman", 0) > 0
    ]

    if not molecules_with_costs:
        return None

    names = [m.get("name") for m in molecules_with_costs]
    costs = [m.get("market", {}).get("total_price_toman", 0)
             for m in molecules_with_costs]

    # Create chart
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(names, costs, color="#667eea", alpha=0.8)

    ax.set_xlabel("Cost (Toman)", fontsize=13)
    ax.set_title("Raw Material Cost Comparison",
                fontsize=16, fontweight='bold', pad=20)
    ax.grid(axis='x', alpha=0.3, linestyle='--')

    # Add value labels
    for bar in bars:
        width = bar.get_width()
        ax.text(width, bar.get_y() + bar.get_height()/2,
               f'{width:,.0f}',
               ha='left', va='center', fontsize=11,
               fontweight='bold', color="#333")

    plt.tight_layout()
    return fig


def get_analyzed_molecules() -> List[str]:
    """Get list of analyzed molecule names.

    Returns:
        List of molecule names
    """
    if not agents.molecule_data:
        return []
    return [m.get("name", "Unknown") for m in agents.molecule_data]


# Build Gradio Interface
with gr.Blocks(
    theme=gr.themes.Soft(
        primary_hue="indigo",
        secondary_hue="purple",
    ),
    title="MoleculeAI - Professional Analysis Platform",
    css="""
        .gradio-container {font-family: 'Inter', sans-serif;}
        .header {background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                 padding: 2rem; border-radius: 10px; color: white;
                 margin-bottom: 2rem; text-align: center;}
        .stat-box {background: #f8f9fa; padding: 1rem; border-radius: 8px;
                   border-left: 4px solid #667eea;}
    """
) as app:

    # Header
    gr.HTML("""
        <div class="header">
            <h1 style="margin:0; font-size: 2.5rem;">🧬 MoleculeAI</h1>
            <p style="margin: 0.5rem 0 0 0; font-size: 1.2rem; opacity: 0.95;">
                Professional Molecular Analysis for Chemical Ingredient Companies
            </p>
        </div>
    """)

    # Main tabs
    with gr.Tabs() as tabs:

        # Tab 1: Analysis
        with gr.Tab("🔬 Molecule Analysis"):
            gr.Markdown("### Upload and Analyze Molecules")
            gr.Markdown("Upload a text file with molecule names (one per line) for comprehensive analysis")

            with gr.Row():
                with gr.Column(scale=1):
                    file_input = gr.File(
                        label="📁 Upload Molecule File",
                        file_types=[".txt"],
                        type="binary"
                    )
                    analyze_btn = gr.Button(
                        "🚀 Analyze Molecules",
                        variant="primary",
                        size="lg"
                    )
                    clear_btn = gr.Button("🗑️ Clear Results")

                with gr.Column(scale=2):
                    status_output = gr.Textbox(
                        label="📊 Analysis Status",
                        lines=8,
                        interactive=False
                    )

            gr.Markdown("### Results")

            with gr.Tabs():
                with gr.Tab("📋 Summary Table"):
                    summary_table = gr.Dataframe(
                        label="Molecule Summary",
                        interactive=False
                    )

                with gr.Tab("📄 Detailed Results"):
                    detailed_output = gr.Markdown()

        # Tab 2: AI Chat
        with gr.Tab("💬 AI Assistant"):
            gr.Markdown("### Ask Questions About Your Molecules")
            gr.Markdown("The AI assistant has full context of all analyzed molecules and can answer questions about synthesis, costs, feasibility, and applications.")

            chatbot = gr.Chatbot(
                label="AI Conversation",
                height=500,
                bubble_full_width=False,
                avatar_images=(None, "🤖")
            )

            with gr.Row():
                chat_input = gr.Textbox(
                    label="Your Question",
                    placeholder="Ask about molecules, synthesis routes, costs, or applications...",
                    scale=4
                )
                chat_btn = gr.Button("Send", variant="primary", scale=1)

            gr.Markdown("**Example questions:**")
            gr.Markdown("""
            - "Which molecule is most cost-effective to produce?"
            - "Compare the synthesis complexity of these molecules"
            - "What equipment do I need for the most complex synthesis?"
            - "Which molecules are suitable for cosmetics?"
            """)

        # Tab 3: Market Analytics
        with gr.Tab("📈 Market Insights"):
            gr.Markdown("### Market Analysis & Trends")

            with gr.Row():
                molecule_selector = gr.Dropdown(
                    label="Select Molecule for Market Trend",
                    choices=[],
                    interactive=True
                )
                trend_btn = gr.Button("📊 Generate Market Trend", variant="primary")

            market_plot = gr.Plot(label="Market Value Projection (2023-2027)")

            gr.Markdown("### Cost Comparison")
            cost_btn = gr.Button("💰 Show Cost Comparison", variant="secondary")
            cost_plot = gr.Plot(label="Raw Material Cost Comparison")

    # Event handlers
    def analyze_and_update(file):
        status, detailed, df = analyze_molecules(file)
        molecules = get_analyzed_molecules()
        return status, detailed, df, gr.Dropdown(choices=molecules)

    def clear_all():
        agents.molecule_data = []
        return "", "", pd.DataFrame(), [], gr.Dropdown(choices=[])

    analyze_btn.click(
        fn=analyze_and_update,
        inputs=[file_input],
        outputs=[status_output, detailed_output, summary_table, molecule_selector]
    )

    clear_btn.click(
        fn=clear_all,
        outputs=[status_output, detailed_output, summary_table, chatbot, molecule_selector]
    )

    chat_btn.click(
        fn=chat_with_ai,
        inputs=[chat_input, chatbot],
        outputs=[chat_input, chatbot]
    )

    chat_input.submit(
        fn=chat_with_ai,
        inputs=[chat_input, chatbot],
        outputs=[chat_input, chatbot]
    )

    trend_btn.click(
        fn=create_market_chart,
        inputs=[molecule_selector],
        outputs=[market_plot]
    )

    cost_btn.click(
        fn=create_cost_comparison,
        outputs=[cost_plot]
    )

    # Footer
    gr.HTML("""
        <div style="text-align: center; margin-top: 2rem; padding: 1rem;
                    border-top: 1px solid #e0e0e0; color: #666;">
            <p><strong>MoleculeAI</strong> - Powered by Google Gemini AI & RDKit Chemistry</p>
            <p style="font-size: 0.9rem;">Built for chemical ingredient companies</p>
        </div>
    """)


if __name__ == "__main__":
    print("🚀 Starting MoleculeAI...")
    print("📍 The app will open in your browser automatically")
    print("🔗 Access at: http://localhost:7860")

    app.launch(
        server_name="0.0.0.0",
        server_port=7860,
        share=False,
        show_error=True,
        quiet=False
    )
