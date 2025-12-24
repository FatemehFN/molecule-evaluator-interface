"""UI components for the application."""
import streamlit as st
import matplotlib.pyplot as plt
from typing import Dict, List


def render_header():
    """Render the application header."""
    st.markdown("""
    <div class="header-container">
        <h1 class="header-title">🧬 MoleculeAI</h1>
        <p class="header-subtitle">
            AI-Powered Molecular Analysis for Chemical Ingredient Companies
        </p>
    </div>
    """, unsafe_allow_html=True)


def render_molecule_card(molecule: Dict):
    """Render a card for a single molecule.

    Args:
        molecule: Dictionary containing molecule analysis data
    """
    with st.container():
        st.markdown(f"""
        <div class="molecule-card">
            <div class="molecule-name">🧪 {molecule.get('name', 'Unknown')}</div>
        </div>
        """, unsafe_allow_html=True)

        # Create tabs for different sections
        tabs = st.tabs(["📊 Properties", "⚗️ Synthesis", "💰 Market", "🏭 Feasibility", "💊 Applications"])

        # Properties tab
        with tabs[0]:
            render_properties_section(molecule)

        # Synthesis tab
        with tabs[1]:
            render_synthesis_section(molecule)

        # Market tab
        with tabs[2]:
            render_market_section(molecule)

        # Feasibility tab
        with tabs[3]:
            render_feasibility_section(molecule)

        # Applications tab
        with tabs[4]:
            render_applications_section(molecule)

        st.markdown("---")


def render_properties_section(molecule: Dict):
    """Render molecular properties section."""
    props = molecule.get("properties")
    smiles = molecule.get("smiles")

    if not props or not smiles:
        st.error("❌ Chemical properties not available")
        return

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Basic Properties**")
        st.info(f"**SMILES:** `{smiles}`")
        st.metric("Molecular Weight", f"{props.get('Molecular Weight', 'N/A')} g/mol")
        st.metric("Molecular Formula", props.get("Molecular Formula", "N/A"))

    with col2:
        st.markdown("**Drug-like Properties**")
        st.metric("LogP", props.get("LogP", "N/A"))
        st.metric("H-Bond Donors", props.get("H-Bond Donors", "N/A"))
        st.metric("H-Bond Acceptors", props.get("H-Bond Acceptors", "N/A"))

    st.metric("Rotatable Bonds", props.get("Rotatable Bonds", "N/A"))
    st.metric("TPSA", f"{props.get('TPSA', 'N/A')} Ų")


def render_synthesis_section(molecule: Dict):
    """Render synthesis information section."""
    synthesis = molecule.get("synthesis", {})

    if not synthesis:
        st.warning("⚠️ Synthesis information not available")
        return

    # Complexity
    complexity = synthesis.get("complexity", "Unknown")
    if "Simple" in complexity:
        st.success(f"**Complexity:** {complexity}")
    elif "Moderate" in complexity:
        st.info(f"**Complexity:** {complexity}")
    else:
        st.warning(f"**Complexity:** {complexity}")

    # Synthesis route
    st.markdown("**Synthesis Route:**")
    route = synthesis.get("synthetic_route", [])
    if route and route != ["Synthesis route not available"]:
        for step in route:
            st.markdown(f"- {step}")
    else:
        st.error("Synthesis route not available")

    # Raw materials
    st.markdown("**Raw Materials:**")
    raw_materials = synthesis.get("raw_materials", [])
    if raw_materials:
        for material in raw_materials:
            st.markdown(f"- {material}")
    else:
        st.info("No raw materials listed")

    # Key intermediates
    intermediates = synthesis.get("key_intermediates", [])
    if intermediates:
        st.markdown("**Key Intermediates:**")
        for intermediate in intermediates:
            st.markdown(f"- {intermediate}")


def render_market_section(molecule: Dict):
    """Render market information section."""
    market = molecule.get("market", {})

    if not market:
        st.warning("⚠️ Market information not available")
        return

    # Price information
    price_info = market.get("price_info", "")
    total_price = market.get("total_price_toman", 0)

    st.markdown("**Raw Material Pricing (Iranian Market):**")
    if price_info and not price_info.startswith("⚠️"):
        st.text(price_info)
        if total_price > 0:
            st.metric("Estimated Total Cost", f"{total_price:,} Toman")
    else:
        st.info("Price information not available")


def render_feasibility_section(molecule: Dict):
    """Render feasibility assessment section."""
    feasibility = molecule.get("feasibility", {})

    if not feasibility:
        st.warning("⚠️ Feasibility assessment not available")
        return

    assessment = feasibility.get("feasibility_assessment", "")

    if "✅" in assessment or "FEASIBLE" in assessment:
        st.success(assessment)
    elif "⚠️" in assessment or "PARTIALLY" in assessment:
        st.warning(assessment)
    elif "❌" in assessment or "NOT FEASIBLE" in assessment:
        st.error(assessment)
    else:
        st.info(assessment)

    # Show equipment availability
    if feasibility.get("equipment_available"):
        st.info("✓ Equipment list loaded and analyzed")
    else:
        st.warning("⚠️ Equipment list not available for analysis")


def render_applications_section(molecule: Dict):
    """Render applications section."""
    applications = molecule.get("applications", "")

    if applications and not applications.startswith("Could not"):
        st.markdown(applications)
    else:
        st.warning("⚠️ Application information not available")


def render_chat_interface(chat_agent, molecule_data: List[Dict], chat_messages: List[Dict]):
    """Render the chat interface.

    Args:
        chat_agent: ChatAgent instance
        molecule_data: List of molecule analysis data
        chat_messages: List of chat message dictionaries
    """
    st.markdown("### 💬 AI Assistant")
    st.markdown("Ask questions about the analyzed molecules, synthesis routes, or market insights.")

    # Chat container
    chat_container = st.container()

    with chat_container:
        # Display chat history
        for message in chat_messages:
            role = message.get("role", "user")
            content = message.get("content", "")

            if role == "user":
                st.markdown(f"""
                <div class="chat-message user-message">
                    <strong>You:</strong><br>{content}
                </div>
                """, unsafe_allow_html=True)
            else:
                st.markdown(f"""
                <div class="chat-message assistant-message">
                    <strong>AI Assistant:</strong><br>{content}
                </div>
                """, unsafe_allow_html=True)

    # Input area
    col1, col2 = st.columns([5, 1])

    with col1:
        user_input = st.text_input(
            "Your question:",
            key="chat_input",
            placeholder="Ask about molecules, synthesis, pricing, or feasibility..."
        )

    with col2:
        send_button = st.button("Send", type="primary", use_container_width=True)

    if send_button and user_input:
        # Add user message
        chat_messages.append({"role": "user", "content": user_input})

        # Get AI response
        with st.spinner("🤔 AI is thinking..."):
            response = chat_agent.process(chat_messages, molecule_data)

        # Add assistant response
        chat_messages.append({"role": "assistant", "content": response})

        # Rerun to update display
        st.rerun()


def render_market_analysis(molecule_data: List[Dict], market_agent):
    """Render market analysis and trends.

    Args:
        molecule_data: List of molecule analysis data
        market_agent: MarketAgent instance
    """
    st.markdown("### 📈 Market Analysis & Trends")

    if not molecule_data:
        st.info("No molecules analyzed yet")
        return

    # Molecule selector
    molecule_names = [m.get("name", "Unknown") for m in molecule_data]
    selected_molecule = st.selectbox(
        "Select a molecule for market trend analysis:",
        molecule_names
    )

    if st.button("📊 Generate Market Trend", type="primary"):
        with st.spinner(f"Analyzing market trends for {selected_molecule}..."):
            trend_data = market_agent.get_market_trend(selected_molecule)

        if trend_data and trend_data.get("global_market"):
            # Create market trend chart
            fig, ax = plt.subplots(figsize=(10, 5))

            years = sorted(trend_data["global_market"].keys())
            global_vals = [trend_data["global_market"][y] for y in years]
            me_vals = [trend_data["middle_east_market"][y] for y in years]

            ax.plot(years, global_vals, marker='o', linewidth=2.5,
                   label="Global Market", color="#667eea")
            ax.plot(years, me_vals, marker='s', linewidth=2.5,
                   label="Middle East Market", color="#764ba2", linestyle="--")

            ax.set_title(f"Market Value Projection: {selected_molecule}",
                        fontsize=16, fontweight='bold')
            ax.set_xlabel("Year", fontsize=12)
            ax.set_ylabel("Market Value (Million USD)", fontsize=12)
            ax.legend(loc='best', fontsize=11)
            ax.grid(True, alpha=0.3)

            st.pyplot(fig)

            # Display data table
            st.markdown("**Market Data Table:**")
            col1, col2 = st.columns(2)

            with col1:
                st.markdown("**Global Market (Million USD)**")
                for year in years:
                    st.metric(year, f"${global_vals[years.index(year)]:.1f}M")

            with col2:
                st.markdown("**Middle East Market (Million USD)**")
                for year in years:
                    st.metric(year, f"${me_vals[years.index(year)]:.1f}M")
        else:
            st.error("Could not retrieve market trend data")

    # Cost comparison
    st.markdown("---")
    st.markdown("### 💰 Cost Comparison")

    molecules_with_costs = [
        m for m in molecule_data
        if m.get("market", {}).get("total_price_toman", 0) > 0
    ]

    if molecules_with_costs:
        names = [m.get("name") for m in molecules_with_costs]
        costs = [m.get("market", {}).get("total_price_toman", 0) for m in molecules_with_costs]

        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.barh(names, costs, color="#667eea")
        ax.set_xlabel("Cost (Toman)", fontsize=12)
        ax.set_title("Raw Material Cost Comparison", fontsize=16, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)

        # Add value labels
        for bar in bars:
            width = bar.get_width()
            ax.text(width, bar.get_y() + bar.get_height()/2,
                   f'{width:,.0f}',
                   ha='left', va='center', fontsize=10)

        st.pyplot(fig)
    else:
        st.info("No cost data available for comparison")
